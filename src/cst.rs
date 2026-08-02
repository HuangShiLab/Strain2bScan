//! Within-species clustering (StrainScan's Cluster Search Tree step) over 2bRAD tags.
//!
//! StrainScan groups near-identical strains into **clusters** before profiling, because
//! strains within ~95% similarity cannot be told apart from short reads — the cluster is
//! the finest reliable resolution. It clusters genomes by k-mer Jaccard (Dashing) with
//! **single-linkage** at **0.05 distance = 0.95 similarity** (`hclsMap_95`). We do the
//! same with 2bRAD tag sets.
//!
//! Single-linkage at threshold τ is exactly the connected components of the graph whose
//! edges are genome pairs with similarity ≥ τ, which we compute with union-find.
//!
//! ## Marker taxonomy (answers: can strain-specific be told from species-specific?)
//!
//! Within a species, each tag is classified by **how many of the species' genomes /
//! clusters carry it** — this is orthogonal to whether the tag is *species*-specific:
//!
//! | class            | cluster degree | genome degree | use |
//! |------------------|----------------|---------------|-----|
//! | `SpeciesCore`    | all clusters   | many          | detect species, NOT strains |
//! | `SharedPartial`  | >1 (not all)   | several       | weak discrimination |
//! | `ClusterSpecific`| exactly 1      | ≥2 in cluster | **Layer-2 cluster marker** |
//! | `StrainSpecific` | exactly 1      | exactly 1     | **finest marker** |
//!
//! The cluster/strain-specific tags are the discriminating markers Layer-2 needs; they are
//! found by within-species incidence, NOT by reusing Fast2bRAD's species-specific DB
//! (which is dominated by `SpeciesCore` and drops tags that recur in other species).

use std::collections::{HashMap, HashSet};

use crate::db::StrainDb;
use crate::fxhash::{FxHashMap, FxHashSet};
use crate::markers::Marker;
use crate::parallel::par_map;

/// Default StrainScan-equivalent similarity cut (0.95 similarity = 0.05 distance).
pub const DEFAULT_SIMILARITY: f64 = 0.95;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MarkerClass {
    SpeciesCore,
    SharedPartial,
    ClusterSpecific(usize),
    StrainSpecific(usize),
}

/// Within-species cluster structure built from genome tag sets.
#[derive(Debug, Clone)]
pub struct SpeciesCst {
    pub genome_names: Vec<String>,
    /// Per genome: single-copy tag markers (used for clustering, scoring, abundance).
    pub genome_markers: Vec<FxHashSet<Marker>>,
    /// Per genome: full tag set (any copy) — used only to define occurrence-based uniqueness.
    pub genome_full: Vec<FxHashSet<Marker>>,
    /// cluster id -> member genome indices.
    pub clusters: Vec<Vec<usize>>,
    /// genome index -> cluster id.
    pub genome_cluster: Vec<usize>,
}

/// Jaccard similarity between two marker sets.
pub fn jaccard(a: &FxHashSet<Marker>, b: &FxHashSet<Marker>) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    let inter = a.iter().filter(|m| b.contains(m)).count();
    let union = a.len() + b.len() - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

/// Max-containment between two marker sets: |A∩B| / min(|A|,|B|). Unlike Jaccard, this stays near
/// 1 when a low-completeness genome's marker set is (approximately) a subset of a complete
/// relative's — so uneven-completeness genomes cluster together instead of fragmenting into spurious
/// singletons. This is the containment estimator used by Mash-screen / sourmash; exact here for
/// small panels. Note max-containment ≥ Jaccard, so containment clustering merges at least as much.
pub fn max_containment(a: &FxHashSet<Marker>, b: &FxHashSet<Marker>) -> f64 {
    if a.is_empty() || b.is_empty() {
        return if a.is_empty() && b.is_empty() { 1.0 } else { 0.0 };
    }
    let inter = a.iter().filter(|m| b.contains(m)).count();
    inter as f64 / a.len().min(b.len()) as f64
}

/// Exact single-linkage on max-containment ≥ threshold (see [`max_containment`]). O(n²·m).
pub fn single_linkage_containment(
    genome_markers: &[FxHashSet<Marker>],
    threshold: f64,
) -> Vec<Vec<usize>> {
    let n = genome_markers.len();
    let rows: Vec<usize> = (0..n).collect();
    let edges: Vec<(usize, usize)> = par_map(&rows, |&i| {
        let mut local = Vec::new();
        for j in (i + 1)..n {
            if max_containment(&genome_markers[i], &genome_markers[j]) >= threshold {
                local.push((i, j));
            }
        }
        local
    })
    .into_iter()
    .flatten()
    .collect();
    union_find_components(n, &edges)
}

/// Connected components of an undirected graph given by `edges` (union-find), as sorted
/// index groups in deterministic order.
fn union_find_components(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for &(i, j) in edges {
        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
        if ri != rj {
            parent[ri] = rj;
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut clusters: Vec<Vec<usize>> = groups.into_values().collect();
    clusters.sort_by_key(|c| c[0]);
    clusters
}

/// Exact single-linkage = connected components of the τ-similarity graph (full-set Jaccard,
/// parallel edge scan). O(n²·m); used for ≲100 genomes.
pub fn single_linkage(genome_markers: &[FxHashSet<Marker>], similarity: f64) -> Vec<Vec<usize>> {
    let n = genome_markers.len();
    let rows: Vec<usize> = (0..n).collect();
    let edges: Vec<(usize, usize)> = par_map(&rows, |&i| {
        let mut local = Vec::new();
        for j in (i + 1)..n {
            if jaccard(&genome_markers[i], &genome_markers[j]) >= similarity {
                local.push((i, j));
            }
        }
        local
    })
    .into_iter()
    .flatten()
    .collect();
    union_find_components(n, &edges)
}

/// Default bottom-k MinHash sketch size for the clustering distance estimator.
pub const SKETCH_K: usize = 2000;
/// Above this many genomes, cluster with MinHash sketches instead of exact Jaccard.
pub const MINHASH_ABOVE: usize = 96;

#[inline]
fn mix64(mut z: u64) -> u64 {
    // splitmix64 finalizer — spreads canonical tag hashes uniformly for MinHash.
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

/// Bottom-k MinHash sketch: the k smallest mixed marker hashes, sorted ascending, deduped.
pub fn minhash_sketch(markers: &FxHashSet<Marker>, k: usize) -> Vec<u64> {
    let mut v: Vec<u64> = markers.iter().map(|&m| mix64(m)).collect();
    v.sort_unstable();
    v.dedup();
    v.truncate(k);
    v
}

/// Estimate Jaccard from two bottom-k sketches (sorted ascending): of the k smallest values
/// in the union, the fraction present in both.
pub fn minhash_jaccard(a: &[u64], b: &[u64], k: usize) -> f64 {
    let (mut i, mut j, mut seen, mut shared) = (0usize, 0usize, 0usize, 0usize);
    while seen < k && (i < a.len() || j < b.len()) {
        if j >= b.len() || (i < a.len() && a[i] < b[j]) {
            i += 1;
        } else if i >= a.len() || (j < b.len() && b[j] < a[i]) {
            j += 1;
        } else {
            shared += 1;
            i += 1;
            j += 1;
        }
        seen += 1;
    }
    if seen == 0 {
        1.0
    } else {
        shared as f64 / seen as f64
    }
}

/// Scalable single-linkage using MinHash sketches (parallel sketch build + parallel edge
/// scan). O(n·m) to sketch + O(n²·k) to compare, with small k ≪ m.
pub fn single_linkage_minhash(
    genome_markers: &[FxHashSet<Marker>],
    similarity: f64,
    k: usize,
) -> Vec<Vec<usize>> {
    let n = genome_markers.len();
    let sketches: Vec<Vec<u64>> = par_map(genome_markers, |s| minhash_sketch(s, k));
    let rows: Vec<usize> = (0..n).collect();
    let edges: Vec<(usize, usize)> = par_map(&rows, |&i| {
        let mut local = Vec::new();
        for j in (i + 1)..n {
            if minhash_jaccard(&sketches[i], &sketches[j], k) >= similarity {
                local.push((i, j));
            }
        }
        local
    })
    .into_iter()
    .flatten()
    .collect();
    union_find_components(n, &edges)
}

/// Scalable single-linkage on max-containment for large panels: the intersection size is estimated
/// from the MinHash-sketch Jaccard and the exact set sizes (|A∩B| = J·(|A|+|B|)/(1+J)), then
/// containment = |A∩B| / min(|A|,|B|). O(n·m) to sketch + O(n²·k) to compare.
pub fn single_linkage_containment_minhash(
    genome_markers: &[FxHashSet<Marker>],
    threshold: f64,
    k: usize,
) -> Vec<Vec<usize>> {
    let n = genome_markers.len();
    let sketches: Vec<Vec<u64>> = par_map(genome_markers, |s| minhash_sketch(s, k));
    let sizes: Vec<usize> = genome_markers.iter().map(|s| s.len()).collect();
    let rows: Vec<usize> = (0..n).collect();
    let edges: Vec<(usize, usize)> = par_map(&rows, |&i| {
        let mut local = Vec::new();
        for j in (i + 1)..n {
            let jac = minhash_jaccard(&sketches[i], &sketches[j], k);
            let (sa, sb) = (sizes[i] as f64, sizes[j] as f64);
            let inter = if jac > 0.0 {
                jac * (sa + sb) / (1.0 + jac)
            } else {
                0.0
            };
            if inter / sa.min(sb).max(1.0) >= threshold {
                local.push((i, j));
            }
        }
        local
    })
    .into_iter()
    .flatten()
    .collect();
    union_find_components(n, &edges)
}

impl SpeciesCst {
    /// Build the CST for one species from `(name, single-copy markers, full tag set)` per
    /// genome. Clustering uses the single-copy markers; the full sets define occurrence-based
    /// uniqueness in `cluster_db`.
    ///
    /// `containment = true` clusters on max-containment instead of Jaccard (see
    /// [`max_containment`]): the fix for reference panels of uneven completeness, where an
    /// incomplete genome's tags are a subset of a complete relative's and Jaccard would spuriously
    /// split them.
    pub fn build(
        genomes: Vec<(String, Vec<Marker>, Vec<Marker>)>,
        similarity: f64,
        containment: bool,
    ) -> Self {
        let genome_names: Vec<String> = genomes.iter().map(|(n, _, _)| n.clone()).collect();
        let genome_full: Vec<FxHashSet<Marker>> =
            genomes.iter().map(|(_, _, f)| f.iter().copied().collect()).collect();
        let genome_markers: Vec<FxHashSet<Marker>> = genomes
            .into_iter()
            .map(|(_, m, _)| m.into_iter().collect())
            .collect();
        // Exact Jaccard for small panels (deterministic, accurate); MinHash sketches scale
        // to large panels where O(n²·m) exact comparison is too costly. STRAIN2BSCAN_CLUSTER
        // = "minhash" | "exact" forces the method (for validation / benchmarking).
        let use_minhash = match std::env::var("STRAIN2BSCAN_CLUSTER").ok().as_deref() {
            Some("minhash") => true,
            Some("exact") => false,
            _ => genome_markers.len() > MINHASH_ABOVE,
        };
        let clusters = match (containment, use_minhash) {
            (true, true) => {
                single_linkage_containment_minhash(&genome_markers, similarity, SKETCH_K)
            }
            (true, false) => single_linkage_containment(&genome_markers, similarity),
            (false, true) => single_linkage_minhash(&genome_markers, similarity, SKETCH_K),
            (false, false) => single_linkage(&genome_markers, similarity),
        };
        let mut genome_cluster = vec![0usize; genome_names.len()];
        for (cid, members) in clusters.iter().enumerate() {
            for &g in members {
                genome_cluster[g] = cid;
            }
        }
        SpeciesCst {
            genome_names,
            genome_markers,
            genome_full,
            clusters,
            genome_cluster,
        }
    }

    pub fn n_clusters(&self) -> usize {
        self.clusters.len()
    }

    /// genome degree (how many genomes carry the tag) and cluster degree.
    fn degrees(&self, m: Marker) -> (usize, HashSet<usize>) {
        let mut gd = 0;
        let mut cset = HashSet::new();
        for (g, set) in self.genome_markers.iter().enumerate() {
            if set.contains(&m) {
                gd += 1;
                cset.insert(self.genome_cluster[g]);
            }
        }
        (gd, cset)
    }

    /// Classify a single tag within this species.
    pub fn classify(&self, m: Marker) -> MarkerClass {
        let (gd, cset) = self.degrees(m);
        let cd = cset.len();
        if cd == self.n_clusters() && self.n_clusters() > 1 {
            MarkerClass::SpeciesCore
        } else if cd == 1 {
            let cid = *cset.iter().next().unwrap();
            if gd == 1 {
                // find the single genome
                let g = (0..self.genome_markers.len())
                    .find(|&g| self.genome_markers[g].contains(&m))
                    .unwrap();
                MarkerClass::StrainSpecific(g)
            } else {
                MarkerClass::ClusterSpecific(cid)
            }
        } else {
            MarkerClass::SharedPartial
        }
    }

    /// Summary counts of each marker class across all tags in the species.
    pub fn marker_class_summary(&self) -> HashMap<&'static str, usize> {
        let mut all: FxHashSet<Marker> = FxHashSet::default();
        for s in &self.genome_markers {
            all.extend(s.iter().copied());
        }
        let mut out: HashMap<&'static str, usize> = HashMap::new();
        for m in all {
            let key = match self.classify(m) {
                MarkerClass::SpeciesCore => "species_core",
                MarkerClass::SharedPartial => "shared_partial",
                MarkerClass::ClusterSpecific(_) => "cluster_specific",
                MarkerClass::StrainSpecific(_) => "strain_specific",
            };
            *out.entry(key).or_insert(0) += 1;
        }
        out
    }

    /// Build a cluster-resolution DB for Layer-2: each cluster is a unit, its markers are
    /// the union of member genomes' tags. `StrainDb::is_unique` on this DB then means
    /// **cluster-specific** automatically (cluster degree == 1).
    pub fn cluster_db(&self) -> StrainDb {
        let units: Vec<(String, Vec<Marker>)> = self
            .clusters
            .iter()
            .enumerate()
            .map(|(cid, members)| {
                let mut set: FxHashSet<Marker> = FxHashSet::default();
                for &g in members {
                    set.extend(self.genome_markers[g].iter().copied());
                }
                // Clean label "C<id>" for TSV joins; genome membership is reported separately.
                (format!("C{cid}"), set.into_iter().collect())
            })
            .collect();
        let mut db = StrainDb::build(units);
        // Occurrence-based uniqueness: a scored (single-copy) marker is unique iff it occurs
        // — at any copy number — in exactly one cluster's genomes. This guards against the
        // single-copy-filter asymmetry that mislabels a tag as cluster-unique when it is
        // multi-copy (hence filtered) in another cluster yet reachable from that cluster's reads.
        let full_unions: Vec<FxHashSet<Marker>> = self
            .clusters
            .iter()
            .map(|members| {
                let mut s = FxHashSet::default();
                for &g in members {
                    s.extend(self.genome_full[g].iter().copied());
                }
                s
            })
            .collect();
        let mut scored: FxHashSet<Marker> = FxHashSet::default();
        for sm in &db.strain_markers {
            scored.extend(sm.iter().copied());
        }
        db.unique_set = scored
            .into_iter()
            .filter(|m| full_unions.iter().filter(|fu| fu.contains(m)).count() == 1)
            .collect();
        db
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two clusters: {g0,g1} share a lot; {g2,g3} share a lot; the two pairs are distant.
    /// (single-copy markers and full set are identical here — synthetic genomes have no
    /// multi-copy tags, so occurrence-based uniqueness reduces to cluster degree.)
    fn two_cluster_species() -> Vec<(String, Vec<Marker>, Vec<Marker>)> {
        let core: Vec<Marker> = (0..200).collect(); // shared by all 4 (species core)
        let clu_a: Vec<Marker> = (200..240).collect(); // g0,g1 only
        let clu_b: Vec<Marker> = (300..340).collect(); // g2,g3 only
        let g = |name: &str, extra: &[Marker], priv_base: Marker| {
            let mut v = core.clone();
            v.extend_from_slice(extra);
            v.extend((0..3).map(|i| priv_base + i)); // 3 private markers
            (name.to_string(), v.clone(), v)
        };
        vec![
            g("g0", &clu_a, 1000),
            g("g1", &clu_a, 1100),
            g("g2", &clu_b, 2000),
            g("g3", &clu_b, 2100),
        ]
    }

    #[test]
    fn clusters_match_expectation_at_0_95() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        assert_eq!(cst.n_clusters(), 2, "clusters: {:?}", cst.clusters);
        // g0 and g1 land in the same cluster; g2/g3 in the other.
        assert_eq!(cst.genome_cluster[0], cst.genome_cluster[1]);
        assert_eq!(cst.genome_cluster[2], cst.genome_cluster[3]);
        assert_ne!(cst.genome_cluster[0], cst.genome_cluster[2]);
    }

    #[test]
    fn containment_clusters_incomplete_genome_with_complete_twin() {
        // An incomplete assembly whose tags are a 50% SUBSET of its complete twin.
        let full: Vec<Marker> = (0..400).collect();
        let partial: Vec<Marker> = (0..200).collect();
        let genomes = vec![
            ("complete".to_string(), full.clone(), full.clone()),
            ("partial".to_string(), partial.clone(), partial),
        ];
        // Jaccard = 200/400 = 0.5 < 0.95 → the incomplete genome spuriously splits off.
        let j = SpeciesCst::build(genomes.clone(), DEFAULT_SIMILARITY, false);
        assert_eq!(
            j.n_clusters(),
            2,
            "Jaccard should split the incomplete genome: {:?}",
            j.clusters
        );
        // Max-containment = 200/min(200,400) = 1.0 ≥ 0.95 → correctly one cluster.
        let c = SpeciesCst::build(genomes, DEFAULT_SIMILARITY, true);
        assert_eq!(
            c.n_clusters(),
            1,
            "containment should merge incomplete with its twin: {:?}",
            c.clusters
        );
        assert_eq!(c.genome_cluster[0], c.genome_cluster[1]);
    }

    #[test]
    fn marker_classification_is_correct() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        assert_eq!(cst.classify(0), MarkerClass::SpeciesCore); // in all 4 / both clusters
                                                               // cluster-A shared marker (200): present in g0,g1 → one cluster, ≥2 genomes
        assert!(matches!(cst.classify(200), MarkerClass::ClusterSpecific(_)));
        // private marker of g0 (1000): one genome
        assert!(matches!(cst.classify(1000), MarkerClass::StrainSpecific(0)));
        let s = cst.marker_class_summary();
        assert_eq!(s["species_core"], 200);
        assert_eq!(s["cluster_specific"], 80); // 40 (clu_a) + 40 (clu_b)
        assert_eq!(s["strain_specific"], 12); // 4 genomes × 3 private
    }

    /// The hierarchy diagnostic must find the group-specific markers a Cluster Search Tree
    /// would store at each internal node. On the two-cluster fixture the node joining g0+g1 is
    /// core-to-that-pair and absent from g2/g3, so it should hold exactly the 40 cluster-A
    /// markers; the root holds the 200 species-core markers.
    #[test]
    fn hierarchy_stats_find_group_specific_markers() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        let stats = cst.hierarchy_stats();
        assert_eq!(stats.len(), 3, "4 genomes -> 3 internal nodes");

        let pair = |s: &NodeStat| -> Vec<usize> { s.member_genomes.clone() };
        let cherry_a = stats.iter().find(|s| pair(s) == vec![0, 1]).expect("g0+g1 node");
        let cherry_b = stats.iter().find(|s| pair(s) == vec![2, 3]).expect("g2+g3 node");
        assert_eq!(cherry_a.group_specific, 40, "cluster-A shared block");
        assert_eq!(cherry_b.group_specific, 40, "cluster-B shared block");

        let root = stats.iter().max_by_key(|s| s.n_members).unwrap();
        assert_eq!(root.n_members, 4);
        assert_eq!(root.group_specific, 200, "root holds the species core");
    }

    /// A star topology — every genome equidistant, no nested structure — must report zero
    /// group-specific markers at every node but the root. This is the case where a tree carries
    /// no signal, and the diagnostic exists to tell them apart.
    #[test]
    fn hierarchy_stats_report_zero_on_a_star_topology() {
        let core: Vec<Marker> = (0..200).collect();
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> = (0..4)
            .map(|i| {
                let mut v = core.clone();
                v.extend((0..50).map(|k| 1000 + i * 100 + k));
                (format!("g{i}"), v.clone(), v)
            })
            .collect();
        let cst = SpeciesCst::build(genomes, DEFAULT_SIMILARITY, false);
        let stats = cst.hierarchy_stats();
        for s in &stats {
            if s.n_members < 4 {
                assert_eq!(s.group_specific, 0, "no nesting -> no group-specific markers");
            }
        }
    }

    #[test]
    fn cluster_db_marks_cluster_specific_as_unique() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        let db = cst.cluster_db();
        assert_eq!(db.n_strains(), 2);
        // a cluster-A marker is unique (cluster-specific) in the cluster DB...
        assert!(db.is_unique(200));
        // ...but the species-core marker is shared across both clusters.
        assert!(!db.is_unique(0));
    }

    /// A LEAF's marker set is the UNION of its members', not the intersection.
    ///
    /// This has to be pinned rather than left to review, because it has already been broken
    /// once: the carrier-set rewrite of `build_tree` originally derived every node the same
    /// way, which silently reimposes the intersection at the leaves. Nothing in the suite
    /// caught it — reverting the union and re-running left all tests green — so the semantics
    /// were resting entirely on someone noticing.
    ///
    /// It matters because `cluster_db` gives a cluster the union of its members' markers. A
    /// tree scoring the same cluster on the intersection is a *weaker* competitor to the flat
    /// path rather than an extension of it, and for a cluster whose members disagree enough
    /// the intersection empties out and the cluster becomes uncallable at any threshold —
    /// measured on real data, a 5-genome C. acnes cluster went from 115 markers to 0.
    #[test]
    fn leaf_marker_set_is_the_union_of_its_members_not_the_intersection() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        let tree = cst.build_tree();

        // g0 and g1 are one cluster; marker 1000 is private to g0 and held by no other genome.
        // Union keeps it (some member carries it, nobody outside does); intersection drops it
        // (g1 does not carry it).
        let leaf = cst.genome_cluster[0];
        assert_eq!(cst.genome_cluster[1], leaf, "g0 and g1 must share a cluster");
        assert!(
            tree.node_markers[leaf].contains(&1000),
            "leaf lost g0's private marker — the intersection has been reimposed at the leaves"
        );
        assert!(
            tree.node_markers[leaf].contains(&1100),
            "leaf lost g1's private marker"
        );
        // The cluster-shared marker is in the union too, so a passing test cannot be explained
        // by the leaf simply holding everything.
        assert!(tree.node_markers[leaf].contains(&200));
        // A marker belonging to the OTHER cluster must not leak in: the leaf is still
        // "minus everything outside".
        assert!(!tree.node_markers[leaf].contains(&300));
        // Nor the species core, which every genome carries.
        assert!(!tree.node_markers[leaf].contains(&0));
    }

    /// An INTERNAL node keeps the intersection — a marker is diagnostic of a clade only if
    /// every genome beneath it carries the marker. This is the other half of the leaf/internal
    /// split, and pinning only the leaf half would let a rewrite unify them in the other
    /// direction.
    #[test]
    fn internal_node_marker_set_is_the_intersection_minus_everything_outside() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY, false);
        let tree = cst.build_tree();
        let leaf_a = cst.genome_cluster[0];

        // The parent of the g0/g1 leaf spans both clusters here (2 clusters -> the root).
        let parent = tree.parent[leaf_a].expect("leaf must have a parent");
        assert!(
            !tree.node_markers[parent].contains(&1000),
            "internal node took a marker only one descendant carries — intersection violated"
        );
        assert!(
            !tree.node_markers[parent].contains(&200),
            "internal node took a marker only one child cluster carries"
        );
        // Its own genomes are every genome, so "minus everything outside" removes nothing and
        // the species core survives as the clade's diagnostic set.
        assert!(tree.node_markers[parent].contains(&0));
    }
}

// ===== Hierarchy feasibility diagnostic ====================================
//
// Before porting StrainScan's Cluster Search Tree we need to know whether a hierarchy carries
// any signal on 2bRAD markers. StrainScan works on the full k-mer set, where an internal node's
// "group-specific" k-mers number in the tens of thousands; 2bRAD tags are ~76-116x sparser, so
// the same nodes could be empty and the descent would degenerate into enumerating every leaf.
//
// This measures it directly rather than arguing from the mechanism. It builds the agglomerative
// hierarchy above the existing clusters and reports, per internal node, the set StrainScan would
// store there:
//
//     K(v) = (intersection of v's descendant genomes) \ (union of all genomes outside v)
//
// which is exactly "core to this subtree and found nowhere else".

/// One internal node of the candidate hierarchy.
#[derive(Debug, Clone)]
pub struct NodeStat {
    pub node_id: usize,
    pub n_members: usize,
    /// Similarity at which the two children merged (max-linkage Jaccard).
    pub merge_similarity: f64,
    /// |intersection of member genomes' markers|.
    pub core: usize,
    /// |K(v)| — the markers StrainScan would store at this node.
    pub group_specific: usize,
    pub member_genomes: Vec<usize>,
}

impl SpeciesCst {
    /// Build the agglomerative hierarchy over this species' genomes and report, for every
    /// internal node, how many group-specific markers it would carry.
    ///
    /// Merging is by **max-linkage** Jaccard between member genomes, matching
    /// `Build_tree.py::hierarchy`'s single-linkage merge on a similarity matrix.
    pub fn hierarchy_stats(&self) -> Vec<NodeStat> {
        let n = self.genome_markers.len();
        if n < 2 {
            return Vec::new();
        }
        // Active nodes: each holds the genome indices beneath it.
        let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let mut active: Vec<usize> = (0..n).collect();
        let mut next_id = n;
        let mut out = Vec::new();

        // Pairwise genome similarity, computed once.
        let sim = |a: usize, b: usize| jaccard(&self.genome_markers[a], &self.genome_markers[b]);
        let mut cache: FxHashMap<(usize, usize), f64> = FxHashMap::default();
        let mut sim_of = |a: usize, b: usize| -> f64 {
            let k = if a < b { (a, b) } else { (b, a) };
            *cache.entry(k).or_insert_with(|| sim(k.0, k.1))
        };

        while active.len() > 1 {
            // Max-linkage: similarity between two nodes is the max over member pairs.
            let (mut best, mut bi, mut bj) = (f64::NEG_INFINITY, 0usize, 1usize);
            for i in 0..active.len() {
                for j in (i + 1)..active.len() {
                    let (a, b) = (active[i], active[j]);
                    let mut s = f64::NEG_INFINITY;
                    for &x in &members[a] {
                        for &y in &members[b] {
                            s = s.max(sim_of(x, y));
                        }
                    }
                    if s > best {
                        best = s;
                        bi = i;
                        bj = j;
                    }
                }
            }
            let (a, b) = (active[bi], active[bj]);
            let mut merged = members[a].clone();
            merged.extend_from_slice(&members[b]);
            merged.sort_unstable();

            // K(v): core to the subtree, absent from every genome outside it.
            let mut core: FxHashSet<Marker> = self.genome_markers[merged[0]].clone();
            for &g in &merged[1..] {
                core.retain(|m| self.genome_markers[g].contains(m));
            }
            let core_n = core.len();
            for g in 0..n {
                if merged.binary_search(&g).is_err() {
                    core.retain(|m| !self.genome_markers[g].contains(m));
                }
            }
            out.push(NodeStat {
                node_id: next_id,
                n_members: merged.len(),
                merge_similarity: best,
                core: core_n,
                group_specific: core.len(),
                member_genomes: merged.clone(),
            });

            members.push(merged);
            active.retain(|&x| x != a && x != b);
            active.push(next_id);
            next_id += 1;
        }
        out
    }
}

// ===== Cluster Search Tree (StrainScan Layer-1) ============================
//
// A real binary hierarchy above the flat clusters, mirroring `Build_tree.py::hierarchy`.
//
// The flat single-linkage partition already IS the τ-cut of this dendrogram, so the leaf set is
// unchanged and the tree is a strict extension: everything the current algorithm does remains
// available, and the tree adds the internal nodes it currently has no way to use.
//
// Each node stores the marker set StrainScan would test there:
//
//     K(v) = (∩ markers of v's descendant genomes) \ (∪ markers of every genome outside v)
//
// i.e. "core to this subtree and found nowhere else in the species". Leaves get their
// cluster-unique markers; internal nodes get the group-specific markers the unique-only
// algorithm discards.

/// A binary hierarchy over a species' clusters, with per-node marker sets.
#[derive(Debug, Clone)]
pub struct Cst {
    /// Leaf id -> member genome indices. Leaves are `0..n_leaves`.
    pub leaves: Vec<Vec<usize>>,
    /// Node id -> parent (None for the root). Length `2*n_leaves - 1`.
    pub parent: Vec<Option<usize>>,
    /// Node id -> its two children (None for leaves).
    pub children: Vec<Option<(usize, usize)>>,
    /// Node id -> the leaf ids beneath it (a leaf's own id for a leaf).
    pub desc_leaves: Vec<Vec<usize>>,
    /// Node id -> K(v).
    pub node_markers: Vec<FxHashSet<Marker>>,
    /// Node id -> similarity at which its children merged (1.0 for leaves).
    pub merge_similarity: Vec<f64>,
    /// Root node id.
    pub root: usize,
}

impl Cst {
    pub fn n_nodes(&self) -> usize {
        self.parent.len()
    }
    pub fn n_leaves(&self) -> usize {
        self.leaves.len()
    }
    #[inline]
    pub fn is_leaf(&self, v: usize) -> bool {
        v < self.leaves.len()
    }
    /// The sibling of `v`, if it has one.
    pub fn sibling(&self, v: usize) -> Option<usize> {
        let p = self.parent[v]?;
        let (a, b) = self.children[p]?;
        Some(if a == v { b } else { a })
    }
}

impl SpeciesCst {
    /// Build the Cluster Search Tree over this species' clusters.
    ///
    /// Merging is by **max-linkage** Jaccard between member genomes — the similarity of two
    /// nodes is the best similarity across their members — matching StrainScan's single-linkage
    /// merge on a similarity matrix. Exactly `n_leaves - 1` internal nodes are created.
    pub fn build_tree(&self) -> Cst {
        let n_leaves = self.clusters.len();
        let n_nodes = if n_leaves == 0 { 0 } else { 2 * n_leaves - 1 };
        let mut parent: Vec<Option<usize>> = vec![None; n_nodes];
        let mut children: Vec<Option<(usize, usize)>> = vec![None; n_nodes];
        let mut desc_leaves: Vec<Vec<usize>> = Vec::with_capacity(n_nodes);
        let mut merge_similarity: Vec<f64> = vec![1.0; n_nodes];
        for l in 0..n_leaves {
            desc_leaves.push(vec![l]);
        }

        let genomes_of = |ls: &[usize], leaves: &[Vec<usize>]| -> Vec<usize> {
            ls.iter().flat_map(|&l| leaves[l].iter().copied()).collect()
        };

        // Cluster-level similarity, computed once. Under single linkage a merged node's
        // similarity to any other is the max of its two children's, so the matrix updates in
        // O(L) per merge — no need to re-derive max-linkage from the member genomes every round,
        // which is what made the build scale as ~O(n^3.6): 0.01 s at 24 genomes but 0.12 s at 48.
        //
        // Filling it is the dominant cost of the whole build on a large panel: it is every
        // genome pair, n(n-1)/2 of them, each an exact Jaccard over tens of thousands of tags,
        // and the merge loop that follows only touches an L x L array of scalars. Measured on
        // 543 genomes, `cluster` spends ~23 s digesting and ~17 s here and in clustering, while
        // dropping from 419 leaves to 61 moves the total by ~3 s — so the row count is not what
        // costs, the pairwise scan is. Rows are independent, so map over them.
        let rows: Vec<usize> = (0..n_leaves).collect();
        let upper: Vec<Vec<f64>> = par_map(&rows, |&i| {
            (0..n_leaves)
                .map(|j| {
                    if j <= i {
                        return f64::NEG_INFINITY;
                    }
                    let mut best = f64::NEG_INFINITY;
                    for &x in &self.clusters[i] {
                        for &y in &self.clusters[j] {
                            best =
                                best.max(jaccard(&self.genome_markers[x], &self.genome_markers[y]));
                        }
                    }
                    best
                })
                .collect()
        });
        let mut sim = vec![vec![f64::NEG_INFINITY; n_leaves]; n_leaves];
        for i in 0..n_leaves {
            for j in (i + 1)..n_leaves {
                sim[i][j] = upper[i][j];
                sim[j][i] = upper[i][j];
            }
        }
        // A merged node reuses its left child's row.
        let mut row_of: Vec<usize> = (0..n_leaves).collect();
        let mut active: Vec<usize> = (0..n_leaves).collect();
        let mut next_id = n_leaves;
        while active.len() > 1 {
            let (mut best, mut bi, mut bj) = (f64::NEG_INFINITY, 0usize, 1usize);
            for i in 0..active.len() {
                for j in (i + 1)..active.len() {
                    let s = sim[row_of[active[i]]][row_of[active[j]]];
                    if s > best {
                        best = s;
                        bi = i;
                        bj = j;
                    }
                }
            }
            let (a, b) = (active[bi], active[bj]);
            let (ra, rb) = (row_of[a], row_of[b]);
            let mut d = desc_leaves[a].clone();
            d.extend_from_slice(&desc_leaves[b]);
            d.sort_unstable();
            desc_leaves.push(d);
            children[next_id] = Some((a, b));
            parent[a] = Some(next_id);
            parent[b] = Some(next_id);
            merge_similarity[next_id] = best;
            for &k in &active {
                if k == a || k == b {
                    continue;
                }
                let rk = row_of[k];
                let v = sim[ra][rk].max(sim[rb][rk]);
                sim[ra][rk] = v;
                sim[rk][ra] = v;
            }
            row_of.push(ra);
            active.retain(|&x| x != a && x != b);
            active.push(next_id);
            next_id += 1;
        }
        let root = if n_nodes == 0 { 0 } else { n_nodes - 1 };

        // --- K(v) for every node, by carrier-set identity.
        //
        // Both cases below ask a question about a marker's CARRIER SET — the genomes that hold
        // it — so indexing markers by that set once answers every node in a single pass, instead
        // of intersecting members and subtracting every outside genome per node
        // (O(nodes x genomes x markers), an L-fold waste).
        //
        //   internal: K(v) = (∩ members) \ (∪ outside)  ⟺  carriers(m) == genomes under v
        //   leaf:     K(v) = (∪ members) \ (∪ outside)  ⟺  carriers(m) ⊆ members(leaf)
        //
        // The leaf case is the union deliberately, not an oversight: `cluster_db` gives a
        // cluster the union of its members' markers, and a tree scoring the same cluster on the
        // intersection would be a weaker competitor to the flat path rather than an extension of
        // it — for a cluster whose members disagree enough the intersection empties and the
        // cluster becomes uncallable at any threshold. Internal nodes keep the intersection,
        // since a marker is diagnostic of a clade only if every genome under it carries it.
        //
        // Because leaves partition the genomes, "carriers ⊆ members(leaf)" is simply "every
        // carrier maps to the same leaf", which needs no set comparison at all.
        let n_genomes = self.genome_markers.len();
        let mut carriers: FxHashMap<Marker, Vec<usize>> = FxHashMap::default();
        for (g, ms) in self.genome_markers.iter().enumerate().take(n_genomes) {
            for &m in ms {
                carriers.entry(m).or_default().push(g);
            }
        }
        let mut node_markers: Vec<FxHashSet<Marker>> = vec![FxHashSet::default(); n_nodes];
        let mut by_exact: FxHashMap<Vec<usize>, FxHashSet<Marker>> = FxHashMap::default();
        for (m, mut gs) in carriers {
            gs.sort_unstable();
            // Leaf: every carrier in one cluster.
            let l0 = self.genome_cluster[gs[0]];
            if gs.iter().all(|&g| self.genome_cluster[g] == l0) {
                node_markers[l0].insert(m);
            }
            // Internal: carrier set exactly equals some node's genome set.
            by_exact.entry(gs).or_default().insert(m);
        }
        for v in n_leaves..n_nodes {
            let mut gs = genomes_of(&desc_leaves[v], &self.clusters);
            gs.sort_unstable();
            if let Some(ms) = by_exact.get(&gs) {
                node_markers[v] = ms.clone();
            }
        }

        Cst {
            leaves: self.clusters.clone(),
            parent,
            children,
            desc_leaves,
            node_markers,
            merge_similarity,
            root,
        }
    }
}
