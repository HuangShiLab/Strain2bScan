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

/// One-per-process warning for legacy database headers (see `load`).
static LEGACY_HEADER_WARN: std::sync::Once = std::sync::Once::new();

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
    /// Cached unique + quantifiable marker panel for each strain. Built on load / restrict so the
    /// hot profiling loop does not repeatedly filter `strain_markers` through `is_unique` and
    /// `is_quantifiable`.
    pub quant_panels: Vec<Vec<Marker>>,
    /// Cached tree node markers after applying `quant_mask`. `None` when no cross-species mask is
    /// in force; built once in `restrict_to` so tree descent stops allocating per node.
    pub masked_node_markers: Option<Vec<Vec<Marker>>>,
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
        db.compute_quant_panels();
        db
    }

    /// Populate `quant_panels` from the current `strain_markers`, `unique_set`, and `quant_mask`.
    fn compute_quant_panels(&mut self) {
        self.quant_panels = self
            .strain_markers
            .iter()
            .map(|set| {
                let mut panel: Vec<Marker> = set
                    .iter()
                    .copied()
                    .filter(|&m| self.is_unique(m) && self.is_quantifiable(m))
                    .collect();
                panel.sort_unstable();
                panel
            })
            .collect();
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
        self.compute_quant_panels();
        self.compute_masked_node_markers(allowed);
    }

    /// Populate `masked_node_markers` by filtering each CST node's marker set with `allowed`.
    fn compute_masked_node_markers(&mut self, allowed: &FxHashSet<Marker>) {
        self.masked_node_markers = self.tree.as_ref().map(|t| {
            t.node_markers
                .iter()
                .map(|set| {
                    set.iter()
                        .copied()
                        .filter(|m| allowed.contains(m))
                        .collect()
                })
                .collect()
        });
    }

    /// The unique markers of strain `j` — cluster-specific within this species, and (when a
    /// cross-species restriction is in force) not shared with any other species in the panel.
    ///
    /// Returns a sorted slice into the precomputed `quant_panels` cache.
    pub fn unique_markers(&self, j: usize) -> &[Marker] {
        &self.quant_panels[j]
    }

    pub fn unique_marker_count(&self, j: usize) -> usize {
        self.quant_panels[j].len()
    }

    // ===== persistence (simple, line-oriented text) ========================
    // Format:
    //   line 1:            "#strain2bscan-db\t<n_strains>\t<enzyme_csv>\t<markers_per_strain_csv>"
    //   optional sections: "#unique\t<marker_hex,...>" and the "#tree"/"#node"/"#leaf" CST block
    //   next n lines:      "<strain_name>\t<marker_hex,marker_hex,...>"
    // Sparse and compact; production would use a binary/bgzf layout.
    //
    // In k-mer mode (`--marker-source kmer`) the <enzyme_csv> position holds the single token
    // `kmer<K>s<S>` (e.g. `kmer31s100`) instead of enzyme names; `load` does not validate the
    // field against the enzyme registry, so such databases load unchanged and `profile`
    // auto-detects the marker source from the token.
    //
    // The 4th header field declares each strain's marker count (in strain order) so `load` can
    // reject a truncated file: a cut that drops whole trailing strain sections mismatches
    // <n_strains>, and a cut inside the last strain line mismatches that strain's count.
    // Databases written before this field existed (3-field header) still load, with strain-count
    // validation only. Known hole: a cut landing mid-token in the final hex number can leave a
    // shorter-but-valid number with the count unchanged — catching that needs a checksum and is
    // out of scope for this text format.

    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        let counts = self
            .strain_markers
            .iter()
            .map(|s| s.len().to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            w,
            "#strain2bscan-db\t{}\t{}\t{counts}",
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
                writeln!(
                    w,
                    "#node\t{v}\t{par}\t{ca}\t{cb}\t{:.6}\t{joined}",
                    t.merge_similarity[v]
                )?;
            }
            for (l, gs) in t.leaves.iter().enumerate() {
                let joined = gs
                    .iter()
                    .map(|g| g.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
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
        // `BufWriter::drop` swallows I/O errors, so without an explicit flush a full disk or
        // interrupted write would produce a truncated database while `save` reports success.
        w.flush()?;
        Ok(())
    }

    pub fn load(path: &Path) -> std::io::Result<Self> {
        let bad = |msg: String| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("{}: {msg}", path.display()),
            )
        };
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
        let mut declared_strains: Option<usize> = None;
        let mut declared_counts: Option<Vec<usize>> = None;
        // Malformed hex used to be silently dropped (`filter_map(... .ok())`), shrinking marker
        // sets with no error; every token must parse now.
        let hexset = |csv: &str, lineno: usize| -> std::io::Result<FxHashSet<Marker>> {
            csv.split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    Marker::from_str_radix(s, 16)
                        .map_err(|_| bad(format!("line {lineno}: malformed hex marker '{s}'")))
                })
                .collect()
        };
        for (lineno, line) in reader.lines().enumerate() {
            let line = line?;
            let lineno = lineno + 1;
            if line.starts_with('#') {
                if line.starts_with("#strain2bscan-db") {
                    let f: Vec<&str> = line.split('\t').collect();
                    let n: usize = f.get(1).and_then(|s| s.parse().ok()).ok_or_else(|| {
                        bad(format!("line {lineno}: malformed header (strain count)"))
                    })?;
                    declared_strains = Some(n);
                    if let Some(csv) = f.get(2) {
                        enzymes = csv
                            .split(',')
                            .filter(|s| !s.is_empty())
                            .map(String::from)
                            .collect();
                    }
                    match f.get(3) {
                        Some(csv) => {
                            let counts: Vec<usize> = csv
                                .split(',')
                                .filter(|s| !s.is_empty())
                                .map(|s| {
                                    s.parse::<usize>().map_err(|_| {
                                        bad(format!(
                                            "line {lineno}: malformed header (marker count '{s}')"
                                        ))
                                    })
                                })
                                .collect::<std::io::Result<_>>()?;
                            if counts.len() != n {
                                return Err(bad(format!(
                                    "line {lineno}: header declares {n} strains but {} marker counts",
                                    counts.len()
                                )));
                            }
                            declared_counts = Some(counts);
                        }
                        None => {
                            // Legacy 3-field header (written before per-strain marker counts):
                            // still loads, but only the strain count can be validated — a
                            // truncation that drops whole trailing strains is undetectable.
                            // Warn once per process: a multi-profile run loads hundreds of
                            // per-species DBs and would be flooded otherwise.
                            LEGACY_HEADER_WARN.call_once(|| {
                                eprintln!(
                                    "warning: {}: legacy database header without per-strain marker \
                                     counts; truncation of legacy databases is not fully detectable \
                                     (rebuild the database to enable it)",
                                    path.display()
                                );
                            });
                        }
                    }
                } else if line.starts_with("#unique") {
                    if let Some(csv) = line.split('\t').nth(1) {
                        unique_set = hexset(csv, lineno)?;
                    }
                } else if line.starts_with("#tree\t") {
                    // A corrupt tree section must be a hard error: the old `unwrap_or(0)`
                    // fallbacks yielded a silent empty tree, quietly degrading `--layer1 cst`
                    // to the flat path.
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() < 4 {
                        return Err(bad(format!("line {lineno}: malformed #tree header")));
                    }
                    let n_leaves: usize = f[1].parse().map_err(|_| {
                        bad(format!(
                            "line {lineno}: malformed #tree leaf count '{}'",
                            f[1]
                        ))
                    })?;
                    let n_nodes: usize = f[2].parse().map_err(|_| {
                        bad(format!(
                            "line {lineno}: malformed #tree node count '{}'",
                            f[2]
                        ))
                    })?;
                    tree_root = f[3].parse().map_err(|_| {
                        bad(format!("line {lineno}: malformed #tree root '{}'", f[3]))
                    })?;
                    tree_parent = vec![None; n_nodes];
                    tree_children = vec![None; n_nodes];
                    tree_markers = vec![FxHashSet::default(); n_nodes];
                    tree_sim = vec![1.0; n_nodes];
                    tree_leaves = vec![Vec::new(); n_leaves];
                    have_tree = true;
                } else if line.starts_with("#node\t") {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() < 7 {
                        return Err(bad(format!("line {lineno}: malformed #node line")));
                    }
                    let v: usize = f[1].parse().map_err(|_| {
                        bad(format!("line {lineno}: malformed #node index '{}'", f[1]))
                    })?;
                    if v >= tree_parent.len() {
                        return Err(bad(format!("line {lineno}: #node index {v} out of range")));
                    }
                    let num = |i: usize| -> std::io::Result<i64> {
                        f[i].parse().map_err(|_| {
                            bad(format!("line {lineno}: malformed #node field '{}'", f[i]))
                        })
                    };
                    let (par, ca, cb) = (num(2)?, num(3)?, num(4)?);
                    tree_parent[v] = (par >= 0).then_some(par as usize);
                    tree_children[v] = (ca >= 0 && cb >= 0).then_some((ca as usize, cb as usize));
                    tree_sim[v] = f[5].parse().map_err(|_| {
                        bad(format!(
                            "line {lineno}: malformed #node similarity '{}'",
                            f[5]
                        ))
                    })?;
                    tree_markers[v] = hexset(f[6], lineno)?;
                } else if line.starts_with("#leaf\t") {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() < 3 {
                        return Err(bad(format!("line {lineno}: malformed #leaf line")));
                    }
                    let l: usize = f[1].parse().map_err(|_| {
                        bad(format!("line {lineno}: malformed #leaf index '{}'", f[1]))
                    })?;
                    if l >= tree_leaves.len() {
                        return Err(bad(format!("line {lineno}: #leaf index {l} out of range")));
                    }
                    tree_leaves[l] = f[2]
                        .split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| {
                            x.parse::<usize>().map_err(|_| {
                                bad(format!("line {lineno}: malformed #leaf member '{x}'"))
                            })
                        })
                        .collect::<std::io::Result<Vec<_>>>()?;
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
                .map(|s| {
                    Marker::from_str_radix(s, 16)
                        .map_err(|_| bad(format!("line {lineno}: malformed hex marker '{s}'")))
                })
                .collect::<std::io::Result<Vec<_>>>()?;
            strains.push((name, markers));
        }
        let n_declared =
            declared_strains.ok_or_else(|| bad("missing #strain2bscan-db header".to_string()))?;
        if strains.len() != n_declared {
            return Err(bad(format!(
                "corrupt or truncated database: header declares {n_declared} strains but {} were parsed",
                strains.len()
            )));
        }
        if let Some(counts) = &declared_counts {
            for (j, ((name, markers), &want)) in strains.iter().zip(counts).enumerate() {
                if markers.len() != want {
                    return Err(bad(format!(
                        "corrupt or truncated database: strain {j} ('{name}') declares {want} markers but {} were parsed",
                        markers.len()
                    )));
                }
            }
        }
        let mut db = StrainDb::build(strains);
        db.enzymes = enzymes;
        db.unique_set = unique_set;
        db.compute_quant_panels();
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
        assert_eq!(db.unique_markers(0).first(), Some(&10));
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
        assert_eq!(
            db.unique_markers(0).len(),
            2,
            "10 and 99 are unique within the species"
        );

        // Panel-wide species-specific markers: 99 is shared with another species, so it is out.
        let specific: FxHashSet<Marker> = [1, 2, 10, 20].into_iter().collect();
        db.restrict_to(&specific);

        let kept: Vec<Marker> = db.unique_markers(0).to_vec();
        assert_eq!(
            kept,
            vec![10],
            "only the genuinely species-specific marker may be used"
        );
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

    /// A k-mer-mode database stores `kmer<K>s<S>` in the header's enzyme position. It must
    /// round-trip intact — including the 4th-field per-strain marker counts — and parse back
    /// as a k-mer token so `profile` can auto-detect the marker source.
    #[test]
    fn kmer_mode_roundtrip() {
        let mut db = toy();
        db.enzymes = vec![crate::markers::kmer_db_token(15, 1)];
        let path = std::env::temp_dir().join("s2bs_kmer_db.tsv");
        db.save(&path).unwrap();

        let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
        let header = text.lines().next().unwrap();
        let f: Vec<&str> = header.split('\t').collect();
        assert_eq!(f.len(), 4, "header must keep the 4-field form: {header}");
        assert_eq!(f[2], "kmer15s1");
        assert_eq!(f[3], "4,4,4", "per-strain marker counts must be declared");

        let back = StrainDb::load(&path).unwrap();
        assert_eq!(back.enzymes, vec!["kmer15s1".to_string()]);
        assert_eq!(
            crate::markers::parse_kmer_db_token(&back.enzymes[0]),
            Some((15, 1))
        );
        assert_eq!(back.n_strains(), 3);
        assert!(back.is_unique(30));
        let _ = std::fs::remove_file(&path);
    }

    /// Regression: a truncated database used to load with NO warning, silently dropping
    /// markers (2204 -> 1321 in the reported case) while profiling still emitted plausible
    /// numbers. The header's declared strain count + per-strain marker counts must now make
    /// `load` fail.
    #[test]
    fn truncated_db_is_rejected() {
        let db = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
            ("C".into(), vec![1, 2, 3, 30, 40]),
        ]);
        let path = std::env::temp_dir().join("s2bs_truncated_db.tsv");
        db.save(&path).unwrap();
        let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();

        // Case 1: cut real content, not just the trailing newline — the last strain loses
        // its final two markers, so its parsed count mismatches the header.
        let (head, last) = text.trim_end_matches('\n').rsplit_once('\n').unwrap();
        let keep = {
            // drop the last two comma-separated tokens of the last strain line
            let i = last.rfind(',').unwrap();
            let i = last[..i].rfind(',').unwrap();
            &last[..i]
        };
        std::fs::write(&path, format!("{head}\n{keep}\n")).unwrap();
        let err = StrainDb::load(&path).unwrap_err();
        assert!(
            err.to_string()
                .contains("declares 5 markers but 3 were parsed"),
            "unexpected error: {err}"
        );

        // Case 2: a cut that drops a whole trailing strain section mismatches the strain count.
        db.save(&path).unwrap();
        let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
        let (head, _) = text.rsplit_once('\n').unwrap();
        let (head, _) = head.rsplit_once('\n').unwrap();
        std::fs::write(&path, format!("{head}\n")).unwrap();
        let err = StrainDb::load(&path).unwrap_err();
        assert!(
            err.to_string()
                .contains("declares 3 strains but 2 were parsed"),
            "unexpected error: {err}"
        );
        let _ = std::fs::remove_file(&path);
    }

    /// Databases written before the header carried per-strain marker counts (3-field header)
    /// must still load — with strain-count validation only, plus a once-per-process warning.
    #[test]
    fn legacy_header_without_marker_count_still_loads() {
        let db = toy();
        let path = std::env::temp_dir().join("s2bs_legacy_db.tsv");
        db.save(&path).unwrap();
        // Strip the 4th header field to synthesize a pre-counts (legacy) database.
        let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
        let mut lines = text.lines();
        let header = lines.next().unwrap();
        assert_eq!(
            header.split('\t').count(),
            4,
            "save must write the counts field"
        );
        let mut legacy = header.split('\t').take(3).collect::<Vec<_>>().join("\t");
        for l in lines {
            legacy.push('\n');
            legacy.push_str(l);
        }
        legacy.push('\n');
        std::fs::write(&path, legacy).unwrap();

        let back = StrainDb::load(&path).unwrap();
        assert_eq!(back.n_strains(), 3);
        assert!(back.is_unique(20));
        assert_eq!(back.unique_marker_count(0), 1);
        let _ = std::fs::remove_file(&path);
    }

    /// Regression: malformed hex tokens were silently discarded by
    /// `filter_map(|s| Marker::from_str_radix(s, 16).ok())`, shrinking marker sets with no
    /// error. A corrupt token must now fail the load with line context.
    #[test]
    fn malformed_hex_marker_is_rejected() {
        let db = toy();
        let path = std::env::temp_dir().join("s2bs_badhex_db.tsv");
        db.save(&path).unwrap();
        let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
        let mut lines: Vec<String> = text.lines().map(String::from).collect();
        // Corrupt the first marker token of the first strain line (line 2; toy() writes no
        // #unique/#tree sections).
        let (name, csv) = lines[1].split_once('\t').unwrap();
        let mut toks: Vec<&str> = csv.split(',').collect();
        toks[0] = "not_hex";
        lines[1] = format!("{name}\t{}", toks.join(","));
        std::fs::write(&path, lines.join("\n") + "\n").unwrap();

        let err = StrainDb::load(&path).unwrap_err();
        assert!(
            err.to_string()
                .contains("line 2: malformed hex marker 'not_hex'"),
            "unexpected error: {err}"
        );
        let _ = std::fs::remove_file(&path);
    }
}
