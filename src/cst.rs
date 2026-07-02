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
    pub genome_markers: Vec<HashSet<Marker>>,
    /// Per genome: full tag set (any copy) — used only to define occurrence-based uniqueness.
    pub genome_full: Vec<HashSet<Marker>>,
    /// cluster id -> member genome indices.
    pub clusters: Vec<Vec<usize>>,
    /// genome index -> cluster id.
    pub genome_cluster: Vec<usize>,
}

/// Jaccard similarity between two marker sets.
pub fn jaccard(a: &HashSet<Marker>, b: &HashSet<Marker>) -> f64 {
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
pub fn single_linkage(genome_markers: &[HashSet<Marker>], similarity: f64) -> Vec<Vec<usize>> {
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
pub fn minhash_sketch(markers: &HashSet<Marker>, k: usize) -> Vec<u64> {
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
    genome_markers: &[HashSet<Marker>],
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

impl SpeciesCst {
    /// Build the CST for one species from `(name, single-copy markers, full tag set)` per
    /// genome. Clustering uses the single-copy markers; the full sets define occurrence-based
    /// uniqueness in `cluster_db`.
    pub fn build(genomes: Vec<(String, Vec<Marker>, Vec<Marker>)>, similarity: f64) -> Self {
        let genome_names: Vec<String> = genomes.iter().map(|(n, _, _)| n.clone()).collect();
        let genome_full: Vec<HashSet<Marker>> =
            genomes.iter().map(|(_, _, f)| f.iter().copied().collect()).collect();
        let genome_markers: Vec<HashSet<Marker>> = genomes
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
        let clusters = if use_minhash {
            single_linkage_minhash(&genome_markers, similarity, SKETCH_K)
        } else {
            single_linkage(&genome_markers, similarity)
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
        let mut all: HashSet<Marker> = HashSet::new();
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
                let mut set: HashSet<Marker> = HashSet::new();
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
        let full_unions: Vec<HashSet<Marker>> = self
            .clusters
            .iter()
            .map(|members| {
                let mut s = HashSet::new();
                for &g in members {
                    s.extend(self.genome_full[g].iter().copied());
                }
                s
            })
            .collect();
        let mut scored: HashSet<Marker> = HashSet::new();
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
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY);
        assert_eq!(cst.n_clusters(), 2, "clusters: {:?}", cst.clusters);
        // g0 and g1 land in the same cluster; g2/g3 in the other.
        assert_eq!(cst.genome_cluster[0], cst.genome_cluster[1]);
        assert_eq!(cst.genome_cluster[2], cst.genome_cluster[3]);
        assert_ne!(cst.genome_cluster[0], cst.genome_cluster[2]);
    }

    #[test]
    fn marker_classification_is_correct() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY);
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

    #[test]
    fn cluster_db_marks_cluster_specific_as_unique() {
        let cst = SpeciesCst::build(two_cluster_species(), DEFAULT_SIMILARITY);
        let db = cst.cluster_db();
        assert_eq!(db.n_strains(), 2);
        // a cluster-A marker is unique (cluster-specific) in the cluster DB...
        assert!(db.is_unique(200));
        // ...but the species-core marker is shared across both clusters.
        assert!(!db.is_unique(0));
    }
}
