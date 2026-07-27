//! StrainScan-style Layer-2 strain resolution + abundance, on 2bRAD-tag markers.
//!
//! 1. **Presence detection by unique markers.** A cluster/strain is called present iff enough
//!    of its *unique* (cluster-specific) markers are observed. Restricting to single-cluster
//!    markers makes detection immune to the shared-marker cross-talk that breaks a greedy
//!    set-cover on conspecific panels.
//! 2. **Absolute depth from unique markers.** Each detected cluster's depth is the
//!    **zero-inclusive** trimmed mean count over its unique-marker panel — see
//!    [`unique_marker_depth`]. Getting this right matters at both scopes: the previous
//!    median-over-detected estimator compressed the ratio between an abundant and a rare
//!    cluster *within* a species too. Because it is absolute (reads per tag), it additionally
//!    lets callers derive a cross-species composition when they want one.
//! 3. **Depth-adaptive gating.** The singleton filter and the coverage floor are functions of
//!    the estimated depth rather than fixed constants, so low-input / high-host samples are
//!    not silently thresholded away. See [`min_count_for`] and [`detectable_fraction`].
//! 4. **Post-filter.** Drop calls below `min_rel_abundance`; renormalize.
//!
//! A non-negative Elastic Net solver ([`nonneg_elastic_net`], StrainScan's
//! `ElasticNet(positive=True)`) is provided but deliberately **not** on the main path: after
//! CST clustering the design columns are already ≥5% dissimilar, presence is decided upstream,
//! and a dense strain×marker matrix does not scale. It is kept for experimentation.

use crate::db::StrainDb;
use crate::markers::MarkerCounts;

/// Depth at or above which a genuine marker is essentially never observed exactly once, so
/// `count == 1` can safely be attributed to sequencing error (StrainScan's singleton rule).
pub const SINGLETON_SAFE_DEPTH: f64 = 3.0;

/// Reciprocal of the fraction of **non-zero** observations winsorized before averaging, to keep
/// collapsed repeats and contamination from inflating the depth estimate (top 1%).
const TRIM_FRACTION: usize = 100;

#[derive(Debug, Clone)]
pub struct Params {
    /// Min number of a strain's unique markers (passing the singleton policy) to call it
    /// present (tag-unit support — StrainScan's `msn`; recalibrate on your data).
    pub min_support_markers: usize,
    /// Min fraction of a strain's unique markers detected (StrainScan 0.7). An absolute floor —
    /// see the comment in [`profile`] for why this one is deliberately not depth-scaled.
    pub min_coverage: f64,
    /// Min relative abundance to keep a call (StrainScan 0.02). Always applied to the
    /// **within-species** fraction [`profile`] computes, including under `multi-profile`.
    pub min_rel_abundance: f64,
    /// Scale the singleton filter to the estimated depth. Turning this off restores the fixed
    /// `count >= 2` behaviour of earlier versions.
    pub adaptive: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            // Recalibrated for sparse 2bRAD markers (full-k-mer StrainScan uses msn*k≈1240
            // k-mers; tag markers are ~50-100x sparser, so the floor is in *tag* units).
            min_support_markers: 10,
            min_coverage: 0.1,
            min_rel_abundance: 0.02,
            adaptive: true,
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
    /// **Absolute** per-tag depth (reads per unique marker), zero-inclusive. Comparable across
    /// species — this is the quantity cross-species composition is built from.
    pub depth: f64,
    /// Relative abundance, normalized over whatever set the caller passed to
    /// [`normalize_by_depth`] (within one species DB after [`profile`]).
    pub rel_abundance: f64,
}

/// Fraction of a marker panel that is *reachable* at per-tag depth `lambda`.
///
/// Under Poisson(λ) sampling a marker is seen at least once with probability `1 − e^(−λ)`, so
/// at λ = 0.1 only 10% of a panel can be detected no matter how good the method is. Gating
/// breadth against a fixed constant therefore rejects genuinely present low-abundance strains;
/// gating against `min_coverage × detectable_fraction(λ)` asks the answerable question
/// ("did we see what was reachable?").
#[inline]
pub fn detectable_fraction(lambda: f64) -> f64 {
    if lambda <= 0.0 {
        0.0
    } else {
        1.0 - (-lambda).exp()
    }
}

/// Minimum per-marker count for a marker to count as evidence, given estimated depth.
///
/// At high depth, `count == 1` is dominated by sequencing error and is filtered. At low depth
/// the opposite holds: under Poisson(λ) the share of *detected* markers seen exactly once is
/// `λ / (e^λ − 1)` — 78% at λ = 0.5 — so the fixed `count >= 2` rule discards most of the
/// signal precisely where signal is scarce. Admitting singletons there costs little precision
/// because sequencing errors generate essentially random tags, which almost never coincide
/// with a *specific* cluster's unique-marker panel; the `min_support_markers` floor still
/// requires many independent hits on that one panel.
#[inline]
pub fn min_count_for(lambda: f64) -> u32 {
    if lambda >= SINGLETON_SAFE_DEPTH {
        2
    } else {
        1
    }
}

/// One pass of per-cluster statistics over a marker panel.
#[derive(Debug, Clone, Copy, Default)]
struct PanelStats {
    /// Panel size (unique markers, or all markers when the cluster has no unique ones).
    panel: usize,
    /// Markers with count >= 1.
    detected1: usize,
    /// Markers with count >= 2.
    detected2: usize,
    /// Zero-inclusive trimmed mean count — the absolute depth estimate.
    depth: f64,
}

/// Compute [`PanelStats`] over cluster `j`'s **unique** markers.
///
/// A cluster with no unique markers (e.g. one whose tag set is a subset of another cluster's)
/// returns an empty panel and is therefore never called. That is deliberate: its every marker
/// is also carried by a co-present relative, so any "evidence" for it is that relative's reads.
/// Measuring it over the shared set instead would report it at the relative's depth, which is
/// exactly the shared-marker cross-talk this module exists to avoid. [`strain_unique_coverage`]
/// keeps a full-marker-set fallback, but only for *reporting* coverage of a cluster that has
/// already been called.
fn panel_stats(db: &StrainDb, counts: &MarkerCounts, j: usize) -> PanelStats {
    let mut obs: Vec<u32> = db
        .unique_markers(j)
        .map(|m| counts.get(&m).copied().unwrap_or(0))
        .collect();
    let panel = obs.len();
    if panel == 0 {
        return PanelStats::default();
    }
    obs.sort_unstable();
    // `obs` is ascending, so `detected*` are suffix lengths.
    let first_ge = |t: u32| obs.partition_point(|&c| c < t);
    let detected1 = panel - first_ge(1);
    let detected2 = panel - first_ge(2);

    // Depth = mean count over the WHOLE panel (zeros included — they are the evidence that the
    // strain is rare), with the top 1% of *non-zero* observations winsorized down to the 99th
    // percentile so collapsed repeats and contamination cannot inflate it.
    //
    // Winsorizing rather than discarding, and taking the fraction of the non-zero observations
    // rather than of the panel, both matter: trimming `panel/100` entries deletes real signal
    // once fewer than 1% of the panel is detected, driving the estimate to zero precisely for
    // the rare strains this estimator exists to measure correctly.
    let depth = if detected1 == 0 {
        0.0
    } else {
        let cap_idx = panel.saturating_sub(1 + detected1 / TRIM_FRACTION);
        let cap = obs[cap_idx] as u64;
        let sum: u64 = obs.iter().map(|&c| (c as u64).min(cap)).sum();
        sum as f64 / panel as f64
    };

    PanelStats {
        panel,
        detected1,
        detected2,
        depth,
    }
}

/// Robust per-strain absolute depth: the **zero-inclusive** trimmed mean count over the
/// strain's unique-marker panel.
///
/// The previous estimator took the median over *detected* markers only, which is severely
/// biased at low abundance: a strain covered at 0.05× has a handful of markers at count 1, so
/// its median is 1 — the same value a strain at 1× reports. Ratios between abundant and rare
/// strains were compressed toward uniform, flattening the whole composition. Averaging over
/// the full panel (zeros included) makes the estimate proportional to true depth; trimming the
/// top 1% keeps the robustness the median was there to provide.
pub fn unique_marker_depth(db: &StrainDb, counts: &MarkerCounts, j: usize) -> f64 {
    panel_stats(db, counts, j).depth
}

/// Coverage = fraction of a strain's unique markers detected (count >= 1).
pub fn strain_unique_coverage(db: &StrainDb, counts: &MarkerCounts, j: usize) -> f64 {
    let st = panel_stats(db, counts, j);
    if st.panel == 0 {
        0.0
    } else {
        st.detected1 as f64 / st.panel as f64
    }
}

/// Detect present clusters/strains by their **unique** markers only.
///
/// Returns `(cluster_index, supporting_marker_count)`.
pub fn detect_present(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<(usize, f64)> {
    let mut out = Vec::new();
    for j in 0..db.n_strains() {
        let st = panel_stats(db, counts, j);
        let min_count = if p.adaptive { min_count_for(st.depth) } else { 2 };
        let detected = if min_count >= 2 {
            st.detected2
        } else {
            st.detected1
        };
        if detected >= p.min_support_markers {
            out.push((j, detected as f64));
        }
    }
    out
}

/// Set `rel_abundance` from absolute `depth` over the given set of calls.
///
/// Call this over **all** calls from all species to get a cross-species composition; call it
/// over one species' calls for a within-species composition. `depth` is in reads-per-tag and
/// single-copy tags are one per cell, so the resulting fractions are cell fractions.
pub fn normalize_by_depth(calls: &mut [StrainCall]) {
    let sum: f64 = calls.iter().map(|c| c.depth).sum();
    if sum > 0.0 {
        for c in calls.iter_mut() {
            c.rel_abundance = c.depth / sum;
        }
    } else if !calls.is_empty() {
        // No depth evidence at all: fall back to uniform rather than emitting a column of
        // zeros, so the documented "sums to 1.0" contract holds for any input. `profile` never
        // reaches this (a call needs detected markers, which implies depth > 0); it exists for
        // external callers.
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

/// Profile one species DB: detect present clusters, estimate absolute depth, gate, and
/// normalize **within this DB** (`rel_abundance` sums to 1.0 across the returned calls).
///
/// Within-species is the scope Layer-2 answers for: given that this species is present, how is
/// it split across strains/clusters? Species-level abundance belongs to the species layer
/// (Fast2bRAD-M), not here.
///
/// To build a **cross-species** composition, do **not** concatenate these fractions — each
/// species sums to 1.0 independently, so a 10x-more-abundant species looks identical to a rare
/// one. Pool the calls and renormalize on the absolute [`StrainCall::depth`] instead, via
/// [`normalize_by_depth`] (which is what `multi-profile`'s `global_abundance` column reports).
pub fn profile(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<StrainCall> {
    let mut calls: Vec<StrainCall> = Vec::new();
    for j in 0..db.n_strains() {
        let st = panel_stats(db, counts, j);
        if st.panel == 0 {
            continue;
        }
        let min_count = if p.adaptive { min_count_for(st.depth) } else { 2 };
        let support = if min_count >= 2 {
            st.detected2
        } else {
            st.detected1
        };
        if support < p.min_support_markers {
            continue;
        }
        // Coverage (breadth) is gated against an ABSOLUTE floor, deliberately.
        //
        // Scaling this gate by `detectable_fraction(depth)` is tempting but provably inert:
        // depth is estimated from the same panel whose breadth is being tested, and when every
        // observed count is 0 or 1 (any depth below ~0.5) `depth <= coverage` identically, so
        // `min_coverage * (1 - e^-depth) < coverage` always holds and the gate can never fire —
        // including for a cluster whose panel was hit only by scattered error tags. This floor
        // is the precision guard that stops a large, sparsely-hit panel from being called;
        // low-depth recall is bought with the singleton policy ([`min_count_for`]) and the
        // Layer-1 species gate instead. Lower `--min-coverage` to trade precision for recall.
        let coverage = st.detected1 as f64 / st.panel as f64;
        if coverage < p.min_coverage {
            continue;
        }
        calls.push(StrainCall {
            strain_index: j,
            name: db.strain_names[j].clone(),
            support: support as f64,
            coverage,
            depth: st.depth,
            rel_abundance: 0.0,
        });
    }
    normalize_by_depth(&mut calls);
    filter_by_abundance(&mut calls, p.min_rel_abundance);
    calls
}

/// Non-negative Elastic Net via cyclic coordinate descent with residual maintenance.
/// Minimizes ½‖Xw − y‖² + α·l1·n·‖w‖₁ + ½·α·(1−l1)·n·‖w‖²  s.t. w ≥ 0.
///
/// Not used by [`profile`] — see the module docs for why.
pub fn nonneg_elastic_net(
    cols: &[Vec<f64>],
    y: &[f64],
    alpha: f64,
    l1_ratio: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let k = cols.len();
    let n = y.len();
    let mut w = vec![0.0; k];
    if n == 0 || k == 0 {
        return w;
    }
    let mut r = y.to_vec(); // residual = y − Xw (w starts at 0)
    let col_sq: Vec<f64> = cols.iter().map(|c| c.iter().map(|v| v * v).sum()).collect();
    let l1 = alpha * l1_ratio * n as f64;
    let l2 = alpha * (1.0 - l1_ratio) * n as f64;

    for _ in 0..max_iter {
        let mut max_dw = 0.0_f64;
        for j in 0..k {
            if col_sq[j] == 0.0 {
                continue;
            }
            // rho = X_j·r + col_sq_j·w_j
            let mut rho = col_sq[j] * w[j];
            for i in 0..n {
                rho += cols[j][i] * r[i];
            }
            // Non-negative soft-threshold update.
            let num = rho - l1;
            let wj = if num > 0.0 {
                num / (col_sq[j] + l2)
            } else {
                0.0
            };
            let dw = wj - w[j];
            if dw != 0.0 {
                for i in 0..n {
                    r[i] -= dw * cols[j][i];
                }
                w[j] = wj;
                max_dw = max_dw.max(dw.abs());
            }
        }
        if max_dw < tol {
            break;
        }
    }
    w
}

/// Naive baseline that mimics `strainscan-rust`: score every strain on **all** its
/// markers (shared included), accept any whose total exceeds a single global-ish
/// threshold, no unique-marker covering. Used by the demo to show over-calling.
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
    use crate::markers::Marker;

    /// Build a conspecific DB: `core` shared by all strains, plus private markers each.
    fn conspecific_db(n_strains: usize, core: usize, private: usize) -> (StrainDb, Vec<Vec<Marker>>) {
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
        let a0 = calls.iter().find(|c| c.strain_index == 0).unwrap().rel_abundance;
        let a2 = calls.iter().find(|c| c.strain_index == 2).unwrap().rel_abundance;
        assert!((a0 - 0.7).abs() < 0.06, "a0={a0}");
        assert!((a2 - 0.3).abs() < 0.06, "a2={a2}");

        let naive = naive_profile(&db, &counts, 1240.0);
        assert_eq!(naive.len(), 4, "naive should over-call all 4: {naive:?}");
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

    /// The depth estimator must not compress the abundance ratio between an abundant and a
    /// rare cluster. Cluster A is at 20 reads/tag (full breadth); cluster B is at 0.3
    /// reads/tag, so 30% of its panel is seen, each exactly once.
    ///
    /// Truth: A = 20/20.3 = 98.5%, B = 0.3/20.3 = 1.48%.
    /// Median-over-detected (the old estimator) reports A = 20, B = 1 → 95.2% / 4.8%, i.e. it
    /// over-states the rare cluster by >3x. The zero-inclusive mean must recover ~1.5%.
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
            "B abundance {} should be ~1.5%, not the ~4.8% the median estimator gave",
            dbc.rel_abundance
        );
    }

    /// The outlier guard must **winsorize the top 1% of non-zero observations**, not discard
    /// `panel/100` entries outright. Discarding deletes genuine signal once breadth drops below
    /// ~10% and drives depth to exactly 0 at breadth <= 1% — the opposite of this estimator's
    /// purpose, and invisible unless the test uses a sparse panel.
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

    /// Winsorizing must still neutralize a collapsed-repeat outlier.
    #[test]
    fn repeat_outliers_do_not_inflate_depth() {
        let panel: Vec<Marker> = (10_000..11_000).collect();
        let db = StrainDb::build(vec![("A".into(), panel.clone())]);
        let mut counts = MarkerCounts::default();
        for &m in &panel {
            counts.insert(m, 20);
        }
        counts.insert(panel[0], 5_000); // one collapsed repeat
        let got = unique_marker_depth(&db, &counts, 0);
        assert!((got - 20.0).abs() < 1e-9, "depth {got} should stay 20.0");
    }

    /// A cluster whose markers are all shared with another cluster has NO unique markers, so it
    /// must never be called: every read supporting it is equally explained by its relative.
    /// Measuring such a cluster over the shared panel reports it at the relative's depth and
    /// invents a phantom 50% call — the shared-marker cross-talk this module exists to prevent.
    #[test]
    fn cluster_with_no_unique_markers_is_never_called() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (10_000..10_500).collect(); // strict subset of A
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B_subset".into(), b)]);
        assert_eq!(db.unique_marker_count(1), 0, "B must have no unique markers");

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20); // only strain A is actually present
        }
        let p = Params {
            min_rel_abundance: 0.0,
            ..Params::default()
        };
        let calls = profile(&db, &counts, &p);
        let names: Vec<&str> = calls.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, vec!["A"], "phantom subset cluster was called: {calls:?}");
        assert!((calls[0].rel_abundance - 1.0).abs() < 1e-9);
    }

    /// The coverage floor is the precision guard against a large panel hit sparsely by stray
    /// tags. It must stay absolute: scaling it by `detectable_fraction(depth)` is inert, because
    /// depth <= coverage whenever every count is 0 or 1, so the gate can never fire.
    #[test]
    fn sparsely_hit_large_panel_is_rejected_by_the_coverage_floor() {
        let panel: Vec<Marker> = (10_000..60_000).collect(); // 50k markers
        let db = StrainDb::build(vec![("GHOST".into(), panel.clone())]);
        let mut counts = MarkerCounts::default();
        for &m in panel.iter().take(1_500) {
            counts.insert(m, 1); // 3% breadth of scattered singletons
        }
        let calls = profile(&db, &counts, &Params::default());
        assert!(
            calls.is_empty(),
            "3% breadth of stray singletons must not produce a call: {calls:?}"
        );
    }

    /// A cluster at depth 0.3 has essentially no markers at count >= 2, so the fixed
    /// singleton filter makes it undetectable. Adaptive gating must recover it.
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
            adaptive: false,
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

    /// Cross-species normalization: two species DBs at 10x different depth must keep that 10x
    /// once pooled — the regression that made every species sum to 1.0 independently.
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
        // depths 40:20:4:2 → 60.6% : 30.3% : 6.1% : 3.0%
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
    fn nnls_recovers_known_coefficients() {
        let cols = vec![vec![1.0, 0.0, 1.0, 2.0], vec![0.0, 1.0, 1.0, 1.0]];
        let y = vec![2.0, 3.0, 5.0, 7.0];
        let w = nonneg_elastic_net(&cols, &y, 0.0, 0.5, 5000, 1e-10);
        assert!((w[0] - 2.0).abs() < 1e-3 && (w[1] - 3.0).abs() < 1e-3, "w={w:?}");
    }
}
