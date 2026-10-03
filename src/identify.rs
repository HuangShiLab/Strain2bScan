//! StrainScan-style Layer-2 strain resolution + abundance, on 2bRAD-tag markers.
//!
//! 1. **Presence detection by unique markers.** A cluster/strain is called present iff enough
//!    of its *unique* (cluster-specific) markers are observed.
//! 2. **Absolute depth from unique markers.** Each detected cluster's depth is the
//!    **zero-inclusive** trimmed mean count over its unique-marker panel.
//! 3. **Depth-adaptive gating.** The singleton filter and the coverage floor are functions of
//!    the estimated depth.
//! 4. **Post-filter.** Drop calls below `min_rel_abundance`; renormalize.

// Re-export the items the rest of the crate uses through `strain2bscan::identify::*`.
pub use crate::depth::{
    detectable_fraction, min_count_for, strain_unique_coverage, support_count, unique_marker_depth,
    PanelStats,
};
pub use crate::detect::{detect_present, panel_stats};
pub use crate::enet::{
    build_l2_design, l2_abundance, nonneg_elastic_net, pre_scan, subset_candidates, L2Design,
    MIN_SUBSET_SHARE,
};
pub use crate::tree::{
    descend_tree, descend_tree_inner, descend_tree_masked, resolve_layer1, tree_utility, TreeCall,
    TreeUtility, MAX_FALLBACK_CLADE, MIN_NODE_MARKERS,
};

use crate::db::StrainDb;
use crate::depth::{marker_panel_evidence, panel_stats as panel_stats_fn};
use crate::markers::{Marker, MarkerCounts};

/// Which Layer-1 (presence detection) to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer1 {
    /// Let the database decide.
    Auto,
    /// Score each cluster independently on its own cluster-unique markers.
    Unique,
    /// Descend the Cluster Search Tree.
    Cst,
}

/// Which Layer-2 (abundance) to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer2 {
    /// Zero-inclusive winsorized mean depth over each cluster's own unique markers.
    Depth,
    /// StrainScan-style joint fit over the shared-marker design matrix.
    Enet,
}

#[derive(Debug, Clone)]
pub struct Params {
    /// Min number of a strain's unique markers to call it present.
    pub min_support_markers: usize,
    /// Min fraction of a strain's unique markers detected.
    pub min_coverage: f64,
    /// Min relative abundance to keep a call.
    pub min_rel_abundance: f64,
    /// Min ratio between consecutive sorted abundances to cut the trace tail.
    pub trace_gap: f64,
    /// Absolute abundance floor applied after [`filter_by_trace_gap`].
    pub trace_floor: f64,
    /// Admit `count == 1` markers as evidence when depth is low.
    pub adaptive_singleton: bool,
    /// Minimum `coverage / (1 − e^(−depth))`.
    pub min_consistency: f64,
    /// Scale the Layer-1 species-marker floor down toward what is reachable at depth.
    pub adaptive_floor: bool,
    /// Presence detection to use.
    pub layer1: Layer1,
    /// Abundance estimator to use.
    pub layer2: Layer2,
    /// ElasticNet penalty for [`Layer2::Enet`].
    pub enet_alpha: f64,
    pub enet_l1_ratio: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            min_support_markers: 8,
            min_coverage: 0.1,
            min_rel_abundance: 0.0,
            trace_gap: 0.0,
            trace_floor: 0.0,
            min_consistency: 0.5,
            adaptive_singleton: true,
            adaptive_floor: true,
            layer1: Layer1::Auto,
            layer2: Layer2::Depth,
            enet_alpha: 0.0,
            enet_l1_ratio: 0.5,
        }
    }
}

#[derive(Debug, Clone)]
pub struct StrainCall {
    pub strain_index: usize,
    pub name: String,
    /// Number of unique markers supporting the call at the depth-adaptive count threshold.
    pub support: f64,
    /// Fraction of the strain's unique markers detected in the sample (breadth).
    pub coverage: f64,
    /// **Absolute** per-tag depth (reads per unique marker).
    pub depth: f64,
    /// Total single-copy tags this cluster carries.
    pub n_markers: usize,
    /// Relative abundance, normalized over whatever set the caller passed.
    pub rel_abundance: f64,
}

/// Set `rel_abundance` from absolute `depth` over the given set of calls.
pub fn normalize_by_depth(calls: &mut [StrainCall]) {
    let sum: f64 = calls.iter().map(|c| c.depth).sum();
    if sum > 0.0 {
        for c in calls.iter_mut() {
            c.rel_abundance = c.depth / sum;
        }
    } else if !calls.is_empty() {
        let share = 1.0 / calls.len() as f64;
        for c in calls.iter_mut() {
            c.rel_abundance = share;
        }
    }
}

/// Drop calls below `min_rel`, renormalize the survivors, and sort by descending abundance.
pub fn filter_by_abundance(calls: &mut Vec<StrainCall>, min_rel: f64) {
    calls.retain(|c| c.rel_abundance >= min_rel);
    let kept: f64 = calls.iter().map(|c| c.rel_abundance).sum();
    if kept > 0.0 {
        for c in calls.iter_mut() {
            c.rel_abundance /= kept;
        }
    }
    calls.sort_by(|a, b| {
        b.rel_abundance
            .partial_cmp(&a.rel_abundance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

/// Drop the trace tail below the largest abundance gap.
pub fn filter_by_trace_gap(calls: &mut Vec<StrainCall>, min_ratio: f64, floor: f64) {
    if min_ratio > 0.0 && calls.len() >= 2 {
        calls.sort_by(|a, b| {
            b.rel_abundance
                .partial_cmp(&a.rel_abundance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let (mut cut, mut best) = (calls.len(), min_ratio);
        for i in 0..calls.len() - 1 {
            let (hi, lo) = (calls[i].rel_abundance, calls[i + 1].rel_abundance);
            if lo > 0.0 {
                let ratio = hi / lo;
                if ratio >= best {
                    best = ratio;
                    cut = i + 1;
                }
            }
        }
        calls.truncate(cut);
    }
    if floor > 0.0 {
        calls.retain(|c| c.rel_abundance >= floor);
    }
}

/// Profile one species DB: detect present clusters, estimate absolute depth, gate, and
/// normalize **within this DB**.
pub fn profile(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<StrainCall> {
    let mut calls: Vec<StrainCall> = match resolve_layer1(db, p) {
        Layer1::Auto | Layer1::Unique => profile_unique(db, counts, p),
        Layer1::Cst => match &db.tree {
            Some(tree) => descend_tree_masked(tree, counts, p, db.masked_node_markers.as_deref())
                .into_iter()
                .filter(|c| c.desc_leaves.iter().all(|&l| l < db.n_strains()))
                .map(|c| {
                    let name = c
                        .desc_leaves
                        .iter()
                        .map(|&l| db.strain_names[l].as_str())
                        .collect::<Vec<_>>()
                        .join("|");
                    let first = c.desc_leaves[0];
                    let n_markers = c
                        .desc_leaves
                        .iter()
                        .map(|&l| db.strain_markers[l].len())
                        .max()
                        .unwrap_or(0);
                    StrainCall {
                        strain_index: first,
                        name,
                        support: c.detected as f64,
                        coverage: c.coverage,
                        depth: c.depth,
                        n_markers,
                        rel_abundance: 0.0,
                    }
                })
                .collect(),
            None => profile_unique(db, counts, p),
        },
    };

    if p.layer2 == Layer2::Enet {
        let called: Vec<usize> = calls.iter().map(|c| c.strain_index).collect();
        let extra = subset_candidates(db, counts, &called, p);
        let idx: Vec<usize> = called.iter().chain(extra.iter()).copied().collect();
        if idx.len() > 1 {
            let design = build_l2_design(db, &idx, counts);
            let mut selected = pre_scan(&design, 15, p.min_support_markers);
            for col in called.len()..idx.len() {
                if !selected.contains(&col) {
                    selected.push(col);
                }
            }
            if !selected.is_empty() {
                let w = l2_abundance(&design, &selected, p.enet_alpha, p.enet_l1_ratio);
                let total: f64 = w.iter().sum();
                let fitted: crate::fxhash::FxHashMap<usize, f64> = selected
                    .iter()
                    .zip(w.iter())
                    .map(|(&col, &depth)| (design.clusters[col], depth))
                    .collect();
                calls.retain(|c| fitted.contains_key(&c.strain_index));
                for c in calls.iter_mut() {
                    c.depth = fitted[&c.strain_index];
                }
                for &j in &extra {
                    let Some(&depth) = fitted.get(&j) else {
                        continue;
                    };
                    if total <= 0.0 || depth / total < MIN_SUBSET_SHARE {
                        continue;
                    }
                    let ms: Vec<Marker> = db.strain_markers[j].iter().copied().collect();
                    let ev = marker_panel_evidence(&ms, counts);
                    if ev.panel == 0 {
                        continue;
                    }
                    let support = support_count(&ev, p.adaptive_singleton);
                    calls.push(StrainCall {
                        strain_index: j,
                        name: db.strain_names[j].clone(),
                        support: support as f64,
                        coverage: ev.detected1 as f64 / ev.panel as f64,
                        depth,
                        n_markers: db.strain_markers[j].len(),
                        rel_abundance: 0.0,
                    });
                }
            }
        }
    }

    normalize_by_depth(&mut calls);
    filter_by_trace_gap(&mut calls, p.trace_gap, p.trace_floor);
    filter_by_abundance(&mut calls, p.min_rel_abundance);
    calls
}

/// The flat Layer-1: score each cluster independently on its own unique markers.
pub fn profile_unique(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<StrainCall> {
    let mut calls: Vec<StrainCall> = Vec::new();
    for j in 0..db.n_strains() {
        let st = panel_stats_fn(db, counts, j);
        if st.panel == 0 {
            continue;
        }
        let support = support_count(&st, p.adaptive_singleton);
        if support < p.min_support_markers {
            continue;
        }
        let coverage = st.detected1 as f64 / st.panel as f64;
        if coverage < p.min_coverage {
            continue;
        }
        let expected_breadth = detectable_fraction(st.depth);
        if expected_breadth > 0.0 && coverage / expected_breadth < p.min_consistency {
            continue;
        }
        calls.push(StrainCall {
            strain_index: j,
            name: db.strain_names[j].clone(),
            support: support as f64,
            coverage,
            depth: st.depth,
            n_markers: db.strain_markers[j].len(),
            rel_abundance: 0.0,
        });
    }
    calls
}

/// Naive baseline that mimics `strainscan-rust`: score every strain on **all** its markers.
pub fn naive_profile(db: &StrainDb, counts: &MarkerCounts, min_score: f64) -> Vec<usize> {
    let mut out = Vec::new();
    for j in 0..db.n_strains() {
        let score: f64 = db.strain_markers[j]
            .iter()
            .map(|&m| counts.get(&m).copied().unwrap_or(0) as f64)
            .sum();
        if score >= min_score {
            out.push(j);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a conspecific DB: `core` shared by all strains, plus private markers each.
    fn conspecific_db(
        n_strains: usize,
        core: usize,
        private: usize,
    ) -> (StrainDb, Vec<Vec<Marker>>) {
        let mut strains = Vec::new();
        let mut privates = Vec::new();
        let core_markers: Vec<Marker> = (0..core as Marker).collect();
        for s in 0..n_strains {
            let base = 1_000_000 + (s * private) as Marker;
            let priv_s: Vec<Marker> = (0..private as Marker).map(|i| base + i).collect();
            let mut all = core_markers.clone();
            all.extend_from_slice(&priv_s);
            strains.push((format!("strain{s}"), all));
            privates.push(priv_s);
        }
        (StrainDb::build(strains), privates)
    }

    /// Sample = mixture {strain → abundance} at depth `d`, plus singleton error markers.
    fn synth_sample(db: &StrainDb, mixture: &[(usize, f64)], depth: f64) -> MarkerCounts {
        let mut present: std::collections::HashMap<Marker, f64> = std::collections::HashMap::new();
        for &(j, ab) in mixture {
            for &m in &db.strain_markers[j] {
                *present.entry(m).or_insert(0.0) += ab;
            }
        }
        let mut counts = MarkerCounts::default();
        for (m, frac) in present {
            let c = (depth * frac).round() as u32;
            if c > 0 {
                counts.insert(m, c);
            }
        }
        for e in 0..50u64 {
            counts.insert(9_000_000 + e, 1);
        }
        counts
    }

    #[test]
    fn resolves_conspecific_mixture_where_naive_overcalls() {
        let (db, _priv) = conspecific_db(4, 200, 50);
        let counts = synth_sample(&db, &[(0, 0.7), (2, 0.3)], 30.0);

        let calls = profile(&db, &counts, &Params::default());
        let mut got: Vec<usize> = calls.iter().map(|c| c.strain_index).collect();
        got.sort();
        assert_eq!(got, vec![0, 2], "calls: {calls:?}");
        let a0 = calls
            .iter()
            .find(|c| c.strain_index == 0)
            .unwrap()
            .rel_abundance;
        let a2 = calls
            .iter()
            .find(|c| c.strain_index == 2)
            .unwrap()
            .rel_abundance;
        assert!((a0 - 0.7).abs() < 0.06, "a0={a0}");
        assert!((a2 - 0.3).abs() < 0.06, "a2={a2}");

        let naive = naive_profile(&db, &counts, 1240.0);
        assert_eq!(naive.len(), 4, "naive should over-call all 4: {naive:?}");
    }

    #[test]
    fn abundance_floor_does_not_delete_correctly_estimated_rare_clusters() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![
            ("dominant".into(), a.clone()),
            ("rare".into(), b.clone()),
        ]);

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 30);
        }
        for &m in b.iter().take(400) {
            counts.insert(m, 1);
        }

        let calls = profile(&db, &counts, &Params::default());
        let names: Vec<&str> = calls.iter().map(|c| c.name.as_str()).collect();
        assert!(
            names.contains(&"rare"),
            "default params dropped the rare cluster: {calls:?}"
        );
        let rare = calls.iter().find(|c| c.name == "rare").unwrap();
        let expected = 0.4 / 30.4;
        assert!(
            (rare.rel_abundance - expected).abs() < 0.005,
            "rare abundance {} should be ~{expected:.4}",
            rare.rel_abundance
        );

        let strict = Params {
            min_rel_abundance: 0.02,
            ..Params::default()
        };
        let strict_names: Vec<String> = profile(&db, &counts, &strict)
            .iter()
            .map(|c| c.name.clone())
            .collect();
        assert_eq!(strict_names, vec!["dominant".to_string()]);
    }

    #[test]
    fn trace_gap_cuts_between_community_and_trace_without_hurting_staggered_mocks() {
        let mut defined = vec![
            StrainCall {
                strain_index: 0,
                name: "A".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 10.0,
                n_markers: 100,
                rel_abundance: 0.45,
            },
            StrainCall {
                strain_index: 1,
                name: "B".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 8.0,
                n_markers: 100,
                rel_abundance: 0.36,
            },
            StrainCall {
                strain_index: 2,
                name: "C".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 2.0,
                n_markers: 100,
                rel_abundance: 0.09,
            },
            StrainCall {
                strain_index: 3,
                name: "trace1".into(),
                support: 10.0,
                coverage: 0.2,
                depth: 0.01,
                n_markers: 100,
                rel_abundance: 0.0002,
            },
            StrainCall {
                strain_index: 4,
                name: "trace2".into(),
                support: 10.0,
                coverage: 0.2,
                depth: 0.01,
                n_markers: 100,
                rel_abundance: 0.0001,
            },
        ];
        filter_by_trace_gap(&mut defined, 10.0, 1e-4);
        let names: Vec<&str> = defined.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["A", "B", "C"],
            "defined community: gap should drop trace tail"
        );

        let mut staggered = vec![
            StrainCall {
                strain_index: 0,
                name: " abundant".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 10.0,
                n_markers: 100,
                rel_abundance: 0.50,
            },
            StrainCall {
                strain_index: 1,
                name: "mid".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 5.0,
                n_markers: 100,
                rel_abundance: 0.25,
            },
            StrainCall {
                strain_index: 2,
                name: "rare".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 2.5,
                n_markers: 100,
                rel_abundance: 0.125,
            },
            StrainCall {
                strain_index: 3,
                name: "very_rare".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 1.25,
                n_markers: 100,
                rel_abundance: 0.0625,
            },
        ];
        filter_by_trace_gap(&mut staggered, 10.0, 1e-4);
        let names: Vec<&str> = staggered.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(
            names,
            vec![" abundant", "mid", "rare", "very_rare"],
            "staggered mock: no 10x gap, all members kept"
        );

        let mut disabled = vec![
            StrainCall {
                strain_index: 0,
                name: "A".into(),
                support: 100.0,
                coverage: 1.0,
                depth: 10.0,
                n_markers: 100,
                rel_abundance: 0.50,
            },
            StrainCall {
                strain_index: 1,
                name: "trace".into(),
                support: 10.0,
                coverage: 0.2,
                depth: 0.01,
                n_markers: 100,
                rel_abundance: 0.0001,
            },
        ];
        filter_by_trace_gap(&mut disabled, 0.0, 0.0);
        assert_eq!(
            disabled.len(),
            2,
            "disabled gap filter should keep everything"
        );
    }

    #[test]
    fn shadow_clusters_are_rejected_but_genuine_rare_ones_are_kept() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B".into(), b.clone())]);

        let names = |p: &Params, counts: &MarkerCounts| -> Vec<String> {
            profile(&db, counts, p)
                .iter()
                .map(|c| c.name.clone())
                .collect()
        };
        let loose = Params {
            min_rel_abundance: 0.0,
            min_consistency: 0.0,
            ..Params::default()
        };
        let default = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };

        let mut shadow = MarkerCounts::default();
        for &m in &a {
            shadow.insert(m, 20);
        }
        for &m in b.iter().take(300) {
            shadow.insert(m, 20);
        }
        assert_eq!(
            names(&loose, &shadow),
            vec!["A", "B"],
            "filter off: both called"
        );
        assert_eq!(
            names(&default, &shadow),
            vec!["A"],
            "default must reject the shadow"
        );
        let calls = profile(&db, &shadow, &default);
        assert!((calls[0].rel_abundance - 1.0).abs() < 1e-9);

        let mut genuine = MarkerCounts::default();
        for &m in &a {
            genuine.insert(m, 20);
        }
        for &m in b.iter().take(330) {
            genuine.insert(m, 1);
        }
        let mut got = names(&default, &genuine);
        got.sort();
        assert_eq!(
            got,
            vec!["A".to_string(), "B".to_string()],
            "a genuinely rare cluster with the same breadth must survive"
        );
    }

    #[test]
    fn tree_pooling_recovers_a_leaf_the_flat_algorithm_misses() {
        use crate::cst::SpeciesCst;
        let core: Vec<Marker> = (0..200).collect();
        let ab: Vec<Marker> = (300..400).collect();
        let cd: Vec<Marker> = (400..500).collect();
        let mk = |extra: &[Marker], uniq: std::ops::Range<Marker>| -> Vec<Marker> {
            let mut v = core.clone();
            v.extend_from_slice(extra);
            v.extend(uniq);
            v
        };
        let ga = mk(&ab, 1000..1005);
        let gb = mk(&ab, 1100..1200);
        let gc = mk(&cd, 1200..1300);
        let gd = mk(&cd, 1300..1400);
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> =
            [("A", &ga), ("B", &gb), ("C", &gc), ("D", &gd)]
                .into_iter()
                .map(|(n, g)| (n.to_string(), g.clone(), g.clone()))
                .collect();
        let cst = SpeciesCst::build(genomes, crate::cst::DEFAULT_SIMILARITY, false);
        assert_eq!(cst.n_clusters(), 4, "each genome should be its own cluster");

        let mut counts = MarkerCounts::default();
        for &m in &ga {
            counts.insert(m, 20);
        }

        let db = cst.cluster_db();
        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let flat: Vec<String> = profile(&db, &counts, &p)
            .iter()
            .map(|c| c.name.clone())
            .collect();
        assert!(
            flat.is_empty(),
            "flat algorithm should miss the sparse leaf, got {flat:?}"
        );

        let tree = cst.build_tree();
        let calls = descend_tree(&tree, &counts, &p);
        assert_eq!(
            calls.len(),
            1,
            "tree should call exactly the present leaf: {calls:?}"
        );
        let call = &calls[0];
        assert_eq!(
            tree.leaves[call.node],
            vec![0],
            "the called leaf must be A (genome 0)"
        );
        assert_eq!(
            call.panel, 305,
            "pooled 5 own + 100 A/B-group + 200 root-core"
        );
        assert!(
            call.path.len() == 3,
            "pooled leaf + 2 ancestors, got {:?}",
            call.path
        );
        assert!((call.coverage - 1.0).abs() < 1e-9);
        assert!((call.depth - 20.0).abs() < 0.5, "depth {}", call.depth);
    }

    #[test]
    fn tree_pooling_stops_when_the_sibling_branch_is_entered() {
        use crate::cst::SpeciesCst;
        let core: Vec<Marker> = (0..200).collect();
        let ab: Vec<Marker> = (300..400).collect();
        let cd: Vec<Marker> = (400..500).collect();
        let mk = |extra: &[Marker], uniq: std::ops::Range<Marker>| -> Vec<Marker> {
            let mut v = core.clone();
            v.extend_from_slice(extra);
            v.extend(uniq);
            v
        };
        let ga = mk(&ab, 1000..1100);
        let gb = mk(&ab, 1100..1200);
        let gc = mk(&cd, 1200..1300);
        let gd = mk(&cd, 1300..1400);
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> =
            [("A", &ga), ("B", &gb), ("C", &gc), ("D", &gd)]
                .into_iter()
                .map(|(n, g)| (n.to_string(), g.clone(), g.clone()))
                .collect();
        let cst = SpeciesCst::build(genomes, crate::cst::DEFAULT_SIMILARITY, false);
        let tree = cst.build_tree();

        let mut counts = MarkerCounts::default();
        for &m in ga.iter().chain(gb.iter()) {
            counts.insert(m, 20);
        }
        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let calls = descend_tree(&tree, &counts, &p);
        assert_eq!(calls.len(), 2, "both leaves present: {calls:?}");
        for c in &calls {
            assert_eq!(
                c.path.len(),
                1,
                "sibling entered -> no ancestor pooling, got path {:?}",
                c.path
            );
            assert_eq!(c.panel, 100, "only the leaf's own 100 exclusive markers");
        }
    }

    #[test]
    fn joint_fit_recovers_a_subset_cluster_the_flat_path_cannot_see() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (10_000..10_500).collect();
        let db = StrainDb::build(vec![
            ("A".into(), a.clone()),
            ("B_subset".into(), b.clone()),
        ]);
        assert_eq!(
            db.unique_marker_count(1),
            0,
            "B has no unique markers by construction"
        );

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, if b.contains(&m) { 15 } else { 10 });
        }

        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let flat: Vec<String> = profile(&db, &counts, &p)
            .iter()
            .map(|c| c.name.clone())
            .collect();
        assert_eq!(
            flat,
            vec!["A".to_string()],
            "flat path cannot see the subset cluster"
        );

        let design = build_l2_design(&db, &[0, 1], &counts);
        assert!(
            design.shared_fraction() > 0.4,
            "half the rows should be shared, got {}",
            design.shared_fraction()
        );
        let w = l2_abundance(&design, &[0, 1], 0.0, 0.5);
        assert!((w[0] - 10.0).abs() < 0.05, "w_A = {} should be 10", w[0]);
        assert!((w[1] - 5.0).abs() < 0.05, "w_B = {} should be 5", w[1]);
    }

    #[test]
    fn enet_reaches_the_subset_cluster_through_profile() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (10_000..10_500).collect();
        let db = StrainDb::build(vec![
            ("A".into(), a.clone()),
            ("B_subset".into(), b.clone()),
        ]);
        assert_eq!(
            db.unique_marker_count(1),
            0,
            "B has no unique markers by construction"
        );

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, if b.contains(&m) { 15 } else { 10 });
        }

        let flat = Params {
            min_rel_abundance: 0.0,
            layer1: Layer1::Unique,
            ..Params::default()
        };
        let got: Vec<String> = profile(&db, &counts, &flat)
            .iter()
            .map(|c| c.name.clone())
            .collect();
        assert_eq!(
            got,
            vec!["A".to_string()],
            "flat path still cannot see the subset cluster"
        );

        let enet = Params {
            layer2: Layer2::Enet,
            ..flat
        };
        let calls = profile(&db, &counts, &enet);
        let names: Vec<&str> = calls.iter().map(|c| c.name.as_str()).collect();
        assert!(
            names.contains(&"B_subset"),
            "enet must now reach the subset cluster through profile(), got {names:?}"
        );
        let b_call = calls.iter().find(|c| c.name == "B_subset").unwrap();
        assert!(
            (b_call.rel_abundance - 1.0 / 3.0).abs() < 0.05,
            "B should be ~1/3 of the composition, got {}",
            b_call.rel_abundance
        );
    }

    #[test]
    fn subset_candidates_reject_an_uncontained_cluster() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let mut c: Vec<Marker> = (10_000..10_400).collect();
        c.extend(90_000..90_600);
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("C".into(), c)]);

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 10);
        }
        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        assert!(
            subset_candidates(&db, &counts, &[0], &p).is_empty(),
            "a cluster with unobserved markers of its own is not a subset candidate"
        );
    }

    #[test]
    fn pre_scan_consumes_the_winners_markers() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let dup: Vec<Marker> = (10_000..10_990).collect();
        let far: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![
            ("A".into(), a.clone()),
            ("A_dup".into(), dup),
            ("Far".into(), far.clone()),
        ]);
        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20);
        }

        let design = build_l2_design(&db, &[0, 1, 2], &counts);
        let chosen = pre_scan(&design, 15, 10);
        assert_eq!(
            chosen.first(),
            Some(&0),
            "A explains the most, picked first"
        );
        assert!(
            !chosen.contains(&1),
            "the 99%-duplicate has almost nothing left to explain once A's markers are consumed"
        );
        assert!(!chosen.contains(&2), "the absent cluster explains nothing");
    }

    #[test]
    fn intermediate_strain_resolves_to_the_clade_not_to_nothing() {
        use crate::cst::SpeciesCst;
        let core: Vec<Marker> = (0..200).collect();
        let ab: Vec<Marker> = (300..500).collect();
        let cd: Vec<Marker> = (500..700).collect();
        let mk = |extra: &[Marker], uniq: std::ops::Range<Marker>| -> Vec<Marker> {
            let mut v = core.clone();
            v.extend_from_slice(extra);
            v.extend(uniq);
            v
        };
        let ga = mk(&ab, 1000..1300);
        let gb = mk(&ab, 1300..1600);
        let gc = mk(&cd, 1600..1900);
        let gd = mk(&cd, 1900..2200);
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> =
            [("A", &ga), ("B", &gb), ("C", &gc), ("D", &gd)]
                .into_iter()
                .map(|(n, g)| (n.to_string(), g.clone(), g.clone()))
                .collect();
        let cst = SpeciesCst::build(genomes, crate::cst::DEFAULT_SIMILARITY, false);
        let tree = cst.build_tree();

        let mut counts = MarkerCounts::default();
        for &m in core.iter().chain(ab.iter()) {
            counts.insert(m, 20);
        }
        for m in (1000..1100).chain(1300..1400) {
            counts.insert(m, 20);
        }

        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let calls = descend_tree(&tree, &counts, &p);
        assert_eq!(calls.len(), 1, "one organism, one call: {calls:?}");
        let c = &calls[0];
        assert!(
            !tree.is_leaf(c.node),
            "must resolve to the clade, not to a leaf: node {} is a leaf",
            c.node
        );
        let mut spanned: Vec<usize> = c
            .desc_leaves
            .iter()
            .flat_map(|&l| tree.leaves[l].clone())
            .collect();
        spanned.sort_unstable();
        assert_eq!(
            spanned,
            vec![0, 1],
            "the clade spanned must be exactly A and B"
        );

        let mut db = cst.cluster_db();
        db.tree = Some(tree);
        let named: Vec<String> = profile(
            &db,
            &counts,
            &Params {
                layer1: Layer1::Cst,
                min_rel_abundance: 0.0,
                ..Params::default()
            },
        )
        .iter()
        .map(|c| c.name.clone())
        .collect();
        assert_eq!(named.len(), 1);
        assert!(
            named[0].contains('|'),
            "a clade-level call must name every leaf it spans, got {}",
            named[0]
        );
    }

    #[test]
    fn singleton_errors_do_not_create_calls() {
        let (db, _) = conspecific_db(3, 100, 40);
        let mut counts = MarkerCounts::default();
        for e in 0..100u64 {
            counts.insert(9_000_000 + e, 1);
        }
        assert!(profile(&db, &counts, &Params::default()).is_empty());
    }

    #[test]
    fn depth_estimator_does_not_flatten_rare_clusters() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B".into(), b.clone())]);

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20);
        }
        for &m in b.iter().take(300) {
            counts.insert(m, 1);
        }

        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let calls = profile(&db, &counts, &p);
        assert_eq!(calls.len(), 2, "both clusters must be called: {calls:?}");

        let da = calls.iter().find(|c| c.name == "A").unwrap();
        let dbc = calls.iter().find(|c| c.name == "B").unwrap();
        assert!((da.depth - 20.0).abs() < 1e-9, "A depth {}", da.depth);
        assert!((dbc.depth - 0.3).abs() < 1e-9, "B depth {}", dbc.depth);
        assert!(
            (dbc.rel_abundance - 0.3 / 20.3).abs() < 1e-6,
            "B abundance {} should be ~1.5%",
            dbc.rel_abundance
        );
    }

    #[test]
    fn sparse_panels_are_not_trimmed_into_underestimates() {
        let panel: Vec<Marker> = (10_000..11_000).collect();
        let db = StrainDb::build(vec![("A".into(), panel.clone())]);
        for detected in [5usize, 10, 11, 20, 30, 50, 100, 300] {
            let mut counts = MarkerCounts::default();
            for &m in panel.iter().take(detected) {
                counts.insert(m, 1);
            }
            let got = unique_marker_depth(&db, &counts, 0);
            let want = detected as f64 / 1000.0;
            assert!(
                (got - want).abs() < 1e-9,
                "breadth {detected}/1000: depth {got} should be {want}"
            );
        }
    }

    #[test]
    fn repeat_outliers_do_not_inflate_depth() {
        let panel: Vec<Marker> = (10_000..11_000).collect();
        let db = StrainDb::build(vec![("A".into(), panel.clone())]);
        let mut counts = MarkerCounts::default();
        for &m in &panel {
            counts.insert(m, 20);
        }
        counts.insert(panel[0], 5_000);
        let got = unique_marker_depth(&db, &counts, 0);
        assert!((got - 20.0).abs() < 1e-9, "depth {got} should stay 20.0");
    }

    #[test]
    fn cluster_with_no_unique_markers_is_never_called() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (10_000..10_500).collect();
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B_subset".into(), b)]);
        assert_eq!(
            db.unique_marker_count(1),
            0,
            "B must have no unique markers"
        );

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20);
        }
        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let calls = profile(&db, &counts, &p);
        let names: Vec<&str> = calls.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["A"],
            "phantom subset cluster was called: {calls:?}"
        );
        assert!((calls[0].rel_abundance - 1.0).abs() < 1e-9);
    }

    #[test]
    fn sparsely_hit_large_panel_is_rejected_by_the_coverage_floor() {
        let panel: Vec<Marker> = (10_000..60_000).collect();
        let db = StrainDb::build(vec![("GHOST".into(), panel.clone())]);
        let mut counts = MarkerCounts::default();
        for &m in panel.iter().take(1_500) {
            counts.insert(m, 1);
        }
        let calls = profile(&db, &counts, &Params::default());
        assert!(
            calls.is_empty(),
            "3% breadth of stray singletons must not produce a call: {calls:?}"
        );
    }

    #[test]
    fn adaptive_gating_recovers_low_depth_clusters() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B".into(), b.clone())]);
        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20);
        }
        for &m in b.iter().take(300) {
            counts.insert(m, 1);
        }

        let fixed = Params {
            adaptive_singleton: false,
            adaptive_floor: false,
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let got: Vec<String> = profile(&db, &counts, &fixed)
            .iter()
            .map(|c| c.name.clone())
            .collect();
        assert_eq!(
            got,
            vec!["A".to_string()],
            "fixed gating should miss the low-depth cluster"
        );

        let adaptive = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        assert_eq!(profile(&db, &counts, &adaptive).len(), 2);
    }

    #[test]
    fn pooled_calls_preserve_cross_species_ratio() {
        let mk = |base: Marker| -> Vec<Marker> { (base..base + 500).collect() };
        let db_a = StrainDb::build(vec![("A0".into(), mk(10_000)), ("A1".into(), mk(20_000))]);
        let db_b = StrainDb::build(vec![("B0".into(), mk(30_000)), ("B1".into(), mk(40_000))]);

        let mut counts = MarkerCounts::default();
        for m in 10_000..10_500 {
            counts.insert(m, 40);
        }
        for m in 20_000..20_500 {
            counts.insert(m, 20);
        }
        for m in 30_000..30_500 {
            counts.insert(m, 4);
        }
        for m in 40_000..40_500 {
            counts.insert(m, 2);
        }

        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let mut pooled = profile(&db_a, &counts, &p);
        pooled.extend(profile(&db_b, &counts, &p));
        normalize_by_depth(&mut pooled);

        let get = |n: &str| pooled.iter().find(|c| c.name == n).unwrap().rel_abundance;
        assert!((get("A0") - 40.0 / 66.0).abs() < 1e-6);
        assert!((get("B1") - 2.0 / 66.0).abs() < 1e-6);
        let species_a: f64 = get("A0") + get("A1");
        let species_b: f64 = get("B0") + get("B1");
        assert!(
            (species_a / species_b - 10.0).abs() < 1e-6,
            "species A must stay 10x species B, got {species_a}/{species_b}"
        );
    }

    #[test]
    fn auto_descends_only_when_the_tree_can_change_an_outcome() {
        let u = |rescuable, informative_nodes, fallback_nodes| TreeUtility {
            rescuable,
            informative_nodes,
            fallback_nodes,
        };
        assert!(u(1, 1, 0).worth_descending());
        assert!(!u(1, 0, 0).worth_descending(), "nothing to pool from");
        assert!(!u(0, 1, 0).worth_descending(), "nothing to rescue");
        assert!(u(0, 0, 1).worth_descending());
        assert!(u(0, 5, 2).worth_descending());
        assert!(!u(0, 0, 0).worth_descending());
    }

    #[test]
    fn auto_falls_back_to_unique_without_a_tree() {
        let db = StrainDb::build(vec![
            ("A".into(), (0..500).collect::<Vec<Marker>>()),
            ("B".into(), (500..1000).collect::<Vec<Marker>>()),
        ]);
        assert!(db.tree.is_none());
        assert!(tree_utility(&db, 8).is_none());
        assert_eq!(resolve_layer1(&db, &Params::default()), Layer1::Unique);
    }

    #[test]
    fn explicit_layer1_is_not_second_guessed_by_auto() {
        let db = StrainDb::build(vec![("A".into(), (0..500).collect::<Vec<Marker>>())]);
        for want in [Layer1::Unique, Layer1::Cst] {
            let p = Params {
                layer1: want,
                ..Params::default()
            };
            assert_eq!(resolve_layer1(&db, &p), want);
        }
    }

    #[test]
    fn tree_utility_counts_rescuable_against_the_support_floor() {
        use crate::cst::{SpeciesCst, DEFAULT_SIMILARITY};
        let core: Vec<Marker> = (0..200).collect();
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> = (0..4u64)
            .map(|i| {
                let mut v = core.clone();
                let n = if i == 0 { 3 } else { 100 };
                v.extend(1000 + i * 1000..1000 + i * 1000 + n);
                (format!("g{i}"), v.clone(), v)
            })
            .collect();
        let cst = SpeciesCst::build(genomes, DEFAULT_SIMILARITY, false);
        let mut db = cst.cluster_db();
        db.tree = Some(cst.build_tree());

        assert_eq!(db.n_strains(), 4, "each genome should form its own cluster");
        assert_eq!(
            tree_utility(&db, 2).expect("tree present").rescuable,
            0,
            "every cluster clears a floor of 2"
        );
        assert_eq!(
            tree_utility(&db, 4).unwrap().rescuable,
            1,
            "only g0's cluster (3 unique markers) is below a floor of 4"
        );
        assert_eq!(
            tree_utility(&db, 101).unwrap().rescuable,
            4,
            "every cluster is below a floor of 101"
        );
    }

    #[test]
    fn nnls_recovers_known_coefficients() {
        let cols = vec![vec![1.0, 0.0, 1.0, 2.0], vec![0.0, 1.0, 1.0, 1.0]];
        let y = vec![2.0, 3.0, 5.0, 7.0];
        let w = nonneg_elastic_net(&cols, &y, 0.0, 0.5, 5000, 1e-10);
        assert!(
            (w[0] - 2.0).abs() < 1e-3 && (w[1] - 3.0).abs() < 1e-3,
            "w={w:?}"
        );
    }
}
