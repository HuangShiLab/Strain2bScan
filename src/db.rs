//! Sparse strain × marker database with unique-marker tracking.
//!
//! Unlike `strainscan-rust` (dense `Array2<u8>` serialized to pretty JSON — tens of GB
//! at real scale), this stores, per strain, only the **set of marker hashes it carries**,
//! plus an inverted index `marker -> #strains` so we can flag markers that are unique to
//! a single strain. Unique markers are StrainScan's discriminating signal and, with
//! 2bRAD tags, are exactly the taxonomy-specific tags Fast2bRAD-M already selects in
//! `build_quan_db.rs` (`taxonomies.len() == 1`) — here applied at strain resolution.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::cst::Cst;
use crate::fxhash::{FxHashMap, FxHashSet};
use crate::markers::Marker;

#[derive(Debug, Default, Clone)]
pub struct StrainDb {
    pub strain_names: Vec<String>,
    /// Per strain: the set of markers it carries.
    pub strain_markers: Vec<FxHashSet<Marker>>,
    /// marker -> number of strains carrying it (inverted-index degree).
    pub marker_degree: FxHashMap<Marker, u32>,
    /// Enzyme set used to build this DB (samples must be digested with the same set).
    pub enzymes: Vec<String>,
    /// Truly-unique markers, defined by **full genome occurrence** (a marker absent — at any
    /// copy number — from every other cluster's genomes), set by CST `cluster_db`. Stricter
    /// than `marker_degree == 1` (single-copy membership), which mislabels a tag as unique when
    /// it is multi-copy in another cluster (single-copy-filter asymmetry) and thus reachable
    /// from that cluster's reads. If empty, uniqueness falls back to `marker_degree`.
    pub unique_set: FxHashSet<Marker>,
    /// Optional **cross-species** restriction on which markers may be used for detection and
    /// quantification, applied by `multi-profile` (see [`StrainDb::restrict_to`]).
    ///
    /// `unique_set` / `marker_degree` only know about *this* species: a tag carried by exactly
    /// one cluster here can still occur in a congener's genomes, and then a co-present congener's
    /// reads land on it and inflate this cluster's depth. In a panel with several species of the
    /// same genus (the norm in mock communities and in saliva) that is a systematic abundance
    /// error, not a rare accident. Runtime-only — never serialized, and `None` for a DB loaded on
    /// its own, since a single DB carries no cross-species information.
    pub quant_mask: Option<FxHashSet<Marker>>,
    /// The Cluster Search Tree, when the database was built by `cluster`.
    ///
    /// Layer-1's tree descent tests the *internal* nodes' marker sets, and those cannot be
    /// recovered from `strain_markers`, which holds only the leaves. Persisting the tree is
    /// therefore what makes `--layer1 cst` usable at profile time. Databases written before this
    /// existed have `None` and fall back to the flat path, so old databases stay readable.
    pub tree: Option<Cst>,
}

impl StrainDb {
    /// Build from `(strain_name, markers)` pairs.
    pub fn build(strains: Vec<(String, Vec<Marker>)>) -> Self {
        let mut db = StrainDb::default();
        for (name, markers) in strains {
            let set: FxHashSet<Marker> = markers.into_iter().collect();
            for &m in &set {
                *db.marker_degree.entry(m).or_insert(0) += 1;
            }
            db.strain_names.push(name);
            db.strain_markers.push(set);
        }
        db
    }

    pub fn n_strains(&self) -> usize {
        self.strain_names.len()
    }

    /// Is `marker` unique to a single cluster? Uses the occurrence-based `unique_set` when set
    /// (CST databases), else falls back to single-copy membership degree.
    #[inline]
    pub fn is_unique(&self, marker: Marker) -> bool {
        if self.unique_set.is_empty() {
            self.marker_degree.get(&marker).copied() == Some(1)
        } else {
            self.unique_set.contains(&marker)
        }
    }

    /// May `marker` be used for detection/quantification? True unless a cross-species
    /// restriction is in force and excludes it.
    #[inline]
    pub fn is_quantifiable(&self, marker: Marker) -> bool {
        match &self.quant_mask {
            Some(allowed) => allowed.contains(&marker),
            None => true,
        }
    }

    /// Restrict detection and quantification to `allowed` — the markers that are specific to
    /// this species across the whole panel of databases being profiled together.
    ///
    /// Only the intersection with this DB's own markers is stored, so the mask stays small.
    pub fn restrict_to(&mut self, allowed: &FxHashSet<Marker>) {
        let kept: FxHashSet<Marker> = self
            .marker_degree
            .keys()
            .copied()
            .filter(|m| allowed.contains(m))
            .collect();
        self.quant_mask = Some(kept);
    }

    /// The unique markers of strain `j` — cluster-specific within this species, and (when a
    /// cross-species restriction is in force) not shared with any other species in the panel.
    pub fn unique_markers(&self, j: usize) -> impl Iterator<Item = Marker> + '_ {
        self.strain_markers[j]
            .iter()
            .copied()
            .filter(move |&m| self.is_unique(m) && self.is_quantifiable(m))
    }

    pub fn unique_marker_count(&self, j: usize) -> usize {
        self.unique_markers(j).count()
    }

    // ===== persistence (simple, line-oriented text) ========================
    // Format:
    //   line 1:            "#strain2bscan-db\t<n_strains>"
    //   next n lines:      "<strain_name>\t<marker_hex,marker_hex,...>"
    // Sparse and compact; production would use a binary/bgzf layout.

    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(
            w,
            "#strain2bscan-db\t{}\t{}",
            self.n_strains(),
            self.enzymes.join(",")
        )?;
        if !self.unique_set.is_empty() {
            let joined = self
                .unique_set
                .iter()
                .map(|m| format!("{m:x}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(w, "#unique\t{joined}")?;
        }
        if let Some(t) = &self.tree {
            writeln!(w, "#tree\t{}\t{}\t{}", t.n_leaves(), t.n_nodes(), t.root)?;
            for v in 0..t.n_nodes() {
                let (ca, cb) = match t.children[v] {
                    Some((a, b)) => (a as i64, b as i64),
                    None => (-1, -1),
                };
                let par = t.parent[v].map(|x| x as i64).unwrap_or(-1);
                let joined = t.node_markers[v]
                    .iter()
                    .map(|m| format!("{m:x}"))
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(w, "#node\t{v}\t{par}\t{ca}\t{cb}\t{:.6}\t{joined}", t.merge_similarity[v])?;
            }
            for (l, gs) in t.leaves.iter().enumerate() {
                let joined = gs.iter().map(|g| g.to_string()).collect::<Vec<_>>().join(",");
                writeln!(w, "#leaf\t{l}\t{joined}")?;
            }
        }
        for (name, markers) in self.strain_names.iter().zip(&self.strain_markers) {
            let joined = markers
                .iter()
                .map(|m| format!("{m:x}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(w, "{name}\t{joined}")?;
        }
        Ok(())
    }

    pub fn load(path: &Path) -> std::io::Result<Self> {
        let reader = BufReader::new(File::open(path)?);
        let mut strains = Vec::new();
        let mut enzymes: Vec<String> = Vec::new();
        let mut unique_set: FxHashSet<Marker> = FxHashSet::default();
        let (mut have_tree, mut tree_root) = (false, 0usize);
        let mut tree_parent: Vec<Option<usize>> = Vec::new();
        let mut tree_children: Vec<Option<(usize, usize)>> = Vec::new();
        let mut tree_markers: Vec<FxHashSet<Marker>> = Vec::new();
        let mut tree_sim: Vec<f64> = Vec::new();
        let mut tree_leaves: Vec<Vec<usize>> = Vec::new();
        let hexset = |csv: &str| {
            csv.split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| Marker::from_str_radix(s, 16).ok())
                .collect::<FxHashSet<Marker>>()
        };
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                if line.starts_with("#strain2bscan-db") {
                    if let Some(csv) = line.split('\t').nth(2) {
                        enzymes = csv.split(',').filter(|s| !s.is_empty()).map(String::from).collect();
                    }
                } else if line.starts_with("#unique") {
                    if let Some(csv) = line.split('\t').nth(1) {
                        unique_set = hexset(csv);
                    }
                } else if line.starts_with("#tree\t") {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() >= 4 {
                        let n_nodes: usize = f[2].parse().unwrap_or(0);
                        tree_root = f[3].parse().unwrap_or(0);
                        tree_parent = vec![None; n_nodes];
                        tree_children = vec![None; n_nodes];
                        tree_markers = vec![FxHashSet::default(); n_nodes];
                        tree_sim = vec![1.0; n_nodes];
                        tree_leaves = vec![Vec::new(); f[1].parse().unwrap_or(0)];
                        have_tree = true;
                    }
                } else if line.starts_with("#node\t") {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() >= 7 {
                        if let Ok(v) = f[1].parse::<usize>() {
                            if v < tree_parent.len() {
                                let par: i64 = f[2].parse().unwrap_or(-1);
                                let ca: i64 = f[3].parse().unwrap_or(-1);
                                let cb: i64 = f[4].parse().unwrap_or(-1);
                                tree_parent[v] = (par >= 0).then_some(par as usize);
                                tree_children[v] =
                                    (ca >= 0 && cb >= 0).then_some((ca as usize, cb as usize));
                                tree_sim[v] = f[5].parse().unwrap_or(1.0);
                                tree_markers[v] = hexset(f[6]);
                            }
                        }
                    }
                } else if line.starts_with("#leaf\t") {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() >= 3 {
                        if let Ok(l) = f[1].parse::<usize>() {
                            if l < tree_leaves.len() {
                                tree_leaves[l] = f[2]
                                    .split(',')
                                    .filter(|x| !x.is_empty())
                                    .filter_map(|x| x.parse::<usize>().ok())
                                    .collect();
                            }
                        }
                    }
                }
                continue;
            }
            if line.is_empty() {
                continue;
            }
            let mut it = line.splitn(2, '\t');
            let name = it.next().unwrap_or("").to_string();
            let markers = it
                .next()
                .unwrap_or("")
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| Marker::from_str_radix(s, 16).ok())
                .collect::<Vec<_>>();
            strains.push((name, markers));
        }
        let mut db = StrainDb::build(strains);
        db.enzymes = enzymes;
        db.unique_set = unique_set;
        if have_tree {
            // `desc_leaves` is derivable from the topology, so it is not serialized.
            let n = tree_parent.len();
            let mut desc: Vec<Vec<usize>> = vec![Vec::new(); n];
            for (l, _) in tree_leaves.iter().enumerate() {
                if l < n {
                    desc[l] = vec![l];
                }
            }
            for v in tree_leaves.len()..n {
                if let Some((a, b)) = tree_children[v] {
                    let mut d = desc[a].clone();
                    d.extend_from_slice(&desc[b]);
                    d.sort_unstable();
                    desc[v] = d;
                }
            }
            db.tree = Some(Cst {
                leaves: tree_leaves,
                parent: tree_parent,
                children: tree_children,
                desc_leaves: desc,
                node_markers: tree_markers,
                merge_similarity: tree_sim,
                root: tree_root,
            });
        }
        Ok(db)
    }

    /// Quick DB stats for the `info`/`build` CLI.
    pub fn stats(&self) -> DbStats {
        let total_markers = self.marker_degree.len();
        let unique_total = self.marker_degree.values().filter(|&&d| d == 1).count();
        let avg_markers = if self.n_strains() == 0 {
            0.0
        } else {
            self.strain_markers.iter().map(|s| s.len()).sum::<usize>() as f64
                / self.n_strains() as f64
        };
        DbStats {
            n_strains: self.n_strains(),
            n_markers: total_markers,
            unique_markers: unique_total,
            avg_markers_per_strain: avg_markers,
            unique_fraction: if total_markers == 0 {
                0.0
            } else {
                unique_total as f64 / total_markers as f64
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct DbStats {
    pub n_strains: usize,
    pub n_markers: usize,
    pub unique_markers: usize,
    pub avg_markers_per_strain: f64,
    pub unique_fraction: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> StrainDb {
        // markers 1,2,3 shared "core"; 10/20/30 are private to A/B/C.
        StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
            ("C".into(), vec![1, 2, 3, 30]),
        ])
    }

    #[test]
    fn unique_markers_are_identified() {
        let db = toy();
        assert!(db.is_unique(10) && db.is_unique(20) && db.is_unique(30));
        assert!(!db.is_unique(1));
        assert_eq!(db.unique_marker_count(0), 1);
        assert_eq!(db.unique_markers(0).next(), Some(10));
    }

    /// Cluster-uniqueness is defined within one species, so a tag can be unique to a cluster
    /// here and still occur in a congener's genomes — where a co-present congener's reads land
    /// on it and inflate this cluster's depth. `restrict_to` removes exactly those markers.
    #[test]
    fn restrict_to_drops_markers_shared_with_another_species() {
        // A carries 10 (private) and 99 (also present in another species' DB).
        let db_unrestricted = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 10, 99]),
            ("B".into(), vec![1, 2, 20]),
        ]);
        let mut db = db_unrestricted.clone();
        assert_eq!(db.unique_markers(0).count(), 2, "10 and 99 are unique within the species");

        // Panel-wide species-specific markers: 99 is shared with another species, so it is out.
        let specific: FxHashSet<Marker> = [1, 2, 10, 20].into_iter().collect();
        db.restrict_to(&specific);

        let kept: Vec<Marker> = db.unique_markers(0).collect();
        assert_eq!(kept, vec![10], "only the genuinely species-specific marker may be used");
        assert!(db.is_quantifiable(10) && !db.is_quantifiable(99));
        // Unrestricted DBs (e.g. single-species `profile`) are unaffected.
        assert!(db_unrestricted.is_quantifiable(99));
    }

    /// The tree must survive a save/load round trip, or `--layer1 cst` silently degrades to the
    /// flat path at profile time with no error.
    #[test]
    fn tree_survives_roundtrip() {
        use crate::cst::{SpeciesCst, DEFAULT_SIMILARITY};
        let core: Vec<Marker> = (0..200).collect();
        let mk = |uniq: std::ops::Range<Marker>| -> Vec<Marker> {
            let mut v = core.clone();
            v.extend(uniq);
            v
        };
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> = (0..4u64)
            .map(|i| {
                let g = mk(1000 + i * 100..1100 + i * 100);
                (format!("g{i}"), g.clone(), g)
            })
            .collect();
        let cst = SpeciesCst::build(genomes, DEFAULT_SIMILARITY, false);
        let mut db = cst.cluster_db();
        db.tree = Some(cst.build_tree());
        let before = db.tree.clone().unwrap();

        let path = std::env::temp_dir().join("s2bs_tree_roundtrip.tsv");
        db.save(&path).unwrap();
        let back = StrainDb::load(&path).unwrap();
        let after = back.tree.expect("tree must survive the round trip");

        assert_eq!(after.n_nodes(), before.n_nodes());
        assert_eq!(after.n_leaves(), before.n_leaves());
        assert_eq!(after.root, before.root);
        assert_eq!(after.parent, before.parent);
        assert_eq!(after.children, before.children);
        assert_eq!(after.leaves, before.leaves);
        for v in 0..before.n_nodes() {
            assert_eq!(
                after.node_markers[v], before.node_markers[v],
                "node {v} marker set changed"
            );
        }
        // desc_leaves is derived on load rather than stored; it must still match.
        assert_eq!(after.desc_leaves, before.desc_leaves);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn roundtrip_save_load() {
        let db = toy();
        let dir = std::env::temp_dir();
        let path = dir.join("strain2bscan_test_db.tsv");
        db.save(&path).unwrap();
        let back = StrainDb::load(&path).unwrap();
        assert_eq!(back.n_strains(), 3);
        assert!(back.is_unique(20));
        let _ = std::fs::remove_file(path);
    }
}
