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
//! ## Why not a regression, as StrainScan uses
//!
//! A non-negative Elastic Net solver ([`nonneg_elastic_net`], StrainScan's
//! `ElasticNet(positive=True)`) is provided but deliberately **not** on the main path. The
//! reason is empirical, and it is the opposite of what one would expect.
//!
//! A regression fits over the *full* marker space, shared markers included, so it uses data this
//! module discards — a real statistical-efficiency advantage. The natural assumption is that its
//! L1 penalty also handles shadow clusters (see [`profile`]) by shrinking a weakly-supported
//! component to zero. **It does not.** Run on the shadow scenario — a strain carrying all of
//! cluster A's distinguishing loci and 30% of cluster B's, at 20x — the solver returns:
//!
//! ```text
//!   alpha = 0      w_A = 18.36   w_B = 4.36    (w_B/w_A = 0.238)
//!   alpha = 0.01   w_A = 18.20   w_B = 4.38    (0.241)
//!   alpha = 0.10   w_A = 16.87   w_B = 4.48    (0.266)   <- stronger penalty is WORSE
//! ```
//!
//! The shadow survives at every penalty, and raising alpha makes it worse: shrinking the
//! dominant coefficient leaves residual on the shared core that the minor one absorbs. On the
//! genuine counterpart (B truly present at 0.4x) the same penalty inflates `w_B` from 0.35 to
//! 0.86 — so alpha degrades both cases at once.
//!
//! The structural reason: 300 markers at count 20 *is* strong evidence in a least-squares sense,
//! and no penalty small enough to leave the real strains alone can remove it. More fundamentally,
//! a least-squares objective sees only the mean — "300 markers at 20x" and "1000 markers at 6x"
//! are the same number to it. What separates them is the **breadth** of the evidence, which the
//! objective discards. That is exactly the quantity the depth–breadth consistency test in
//! [`profile`] reads, and why it succeeds where the penalty cannot: on the same scenario it takes
//! the shadow's share of the species from 28% (this module's estimator alone) or ~20% (the
//! regression) to 0%.
//!
//! Abundance accuracy is otherwise a wash — on ATCC MSA-1002 the L1 error to truth is 0.434 here
//! against 0.425 for StrainScan — so the regression's efficiency edge is offset by this module's
//! robustness to outliers (winsorization), its absolute and cross-comparable units, and its
//! immunity to cross-species cross-talk.
//!
//! One case would still favour a regression: a sample strain that is a genuine *mixture* of two
//! references, which should be apportioned rather than called as one or both. If that turns up,
//! the right scope is a **per-species** fit over the species-specific marker space — which does
//! scale (20 clusters × 20 000 markers is ~3 MB), contrary to an earlier note here that judged
//! it by a global matrix.

use crate::db::StrainDb;
use crate::markers::{Marker, MarkerCounts};

/// Depth at or above which a genuine marker is essentially never observed exactly once, so
/// `count == 1` can safely be attributed to sequencing error (StrainScan's singleton rule).
pub const SINGLETON_SAFE_DEPTH: f64 = 3.0;

/// Reciprocal of the fraction of **non-zero** observations winsorized before averaging, to keep
/// collapsed repeats and contamination from inflating the depth estimate (top 1%).
const TRIM_FRACTION: usize = 100;

/// Which Layer-1 (presence detection) to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer1 {
    /// Score each cluster independently on its own cluster-unique markers (the default, and
    /// what every benchmark in this repo was measured with).
    Unique,
    /// Descend the Cluster Search Tree, pooling ancestor markers along the unique path
    /// ([`descend_tree`]). Requires a database built by `cluster` with the tree persisted.
    Cst,
}

/// Which Layer-2 (abundance) to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layer2 {
    /// Zero-inclusive winsorized mean depth over each cluster's own unique markers.
    Depth,
    /// StrainScan-style joint fit: residual pre-scan then non-negative ElasticNet over the
    /// shared-marker design matrix ([`build_l2_design`]).
    Enet,
}

#[derive(Debug, Clone)]
pub struct Params {
    /// Min number of a strain's unique markers (passing the singleton policy) to call it
    /// present (tag-unit support — StrainScan's `msn`; recalibrate on your data).
    pub min_support_markers: usize,
    /// Min fraction of a strain's unique markers detected (StrainScan 0.7). An absolute floor —
    /// see the comment in [`profile`] for why this one is deliberately not depth-scaled.
    pub min_coverage: f64,
    /// Min relative abundance to keep a call, applied to the **within-species** fraction
    /// [`profile`] computes (including under `multi-profile`). **Defaults to 0** — see
    /// [`Params::default`] for why StrainScan's 0.02 is wrong here.
    pub min_rel_abundance: f64,
    /// Admit `count == 1` markers as evidence when the estimated depth is low enough that
    /// genuine markers are mostly singletons (see [`min_count_for`]). Off ⇒ always `count >= 2`.
    ///
    /// Kept **separate** from [`Params::adaptive_floor`] because the two relaxations trade off
    /// differently: admitting singletons buys real sensitivity on sparse 2bRAD panels, but on a
    /// dense multi-enzyme digest against a large species panel it is also the main way an
    /// absent species accumulates spurious evidence.
    pub adaptive_singleton: bool,
    /// Minimum `coverage / (1 − e^(−depth))`. Rejects "shadow" clusters whose observed markers
    /// are far deeper than their breadth allows — see the note in [`profile`]. Set to 0 to
    /// disable. A genuinely present cluster scores ~1 at any depth.
    pub min_consistency: f64,
    /// Scale the Layer-1 species-marker floor down toward what is reachable at the estimated
    /// depth. Off ⇒ the configured floor applies at every depth.
    ///
    /// ⚠️ This relaxation keys on the species' **own** estimated depth, so an *absent* species —
    /// whose depth is near zero — receives the largest relaxation of all. That is backwards, and
    /// it is why large panels lose specificity: see the note on `MIN_FLOOR_FRACTION` in `main`.
    pub adaptive_floor: bool,
    /// Presence detection to use. Defaults to [`Layer1::Unique`].
    pub layer1: Layer1,
    /// Abundance estimator to use. Defaults to [`Layer2::Depth`].
    pub layer2: Layer2,
    /// ElasticNet penalty for [`Layer2::Enet`]. 0 = pure non-negative least squares, which the
    /// shadow experiment in the module docs found strictly better than any positive value.
    pub enet_alpha: f64,
    pub enet_l1_ratio: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            // Recalibrated for sparse 2bRAD markers (full-k-mer StrainScan uses msn*k≈1240
            // k-mers; tag markers are ~50-100x sparser, so the floor is in *tag* units).
            min_support_markers: 10,
            min_coverage: 0.1,
            // 0, NOT StrainScan's 0.02 — and this must stay tied to the depth estimator.
            //
            // A 0.02 floor was survivable only while `depth` came from the median over
            // *detected* markers, which inflates a rare cluster roughly to 1 read/tag no matter
            // how rare it truly is. Removing that bias (see `unique_marker_depth`) made rare
            // clusters report their real, correctly tiny share — and the unchanged 2% floor then
            // deleted them. Measured on a 30x/0.5x mixture: the biased estimator reported the
            // rare cluster at 3.23% (true 1.64%) and it survived; the unbiased one reports 1.57%
            // and 0.02 drops the call entirely. Recall collapsed on exactly the samples full of
            // rare strains — staggered mocks and high-host dilutions.
            //
            // Precision does not depend on this floor: it is carried by `min_support_markers`,
            // `min_coverage` and (in `multi-profile`) the Layer-1 species gate. On the mocks,
            // moving it between 0.02 and 0 changes no false positive, only true positives.
            // Downstream evaluation applies its own presence threshold, so a second hidden one
            // here is redundant. Pass `--min-abundance 0.02` to restore the old behaviour.
            min_rel_abundance: 0.0,
            // 0.5 sits in the middle of a measured gap: synthetic shadows scored <= 0.897 (and
            // real ones far lower), genuine clusters >= 0.949 across depths 0.3x-20x. Left low
            // enough to tolerate the coverage non-uniformity of real 2bRAD, which the synthetic
            // sweep does not model.
            min_consistency: 0.5,
            adaptive_singleton: true,
            adaptive_floor: true,
            layer1: Layer1::Unique,
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
    /// **Absolute** per-tag depth (reads per unique marker), zero-inclusive. Comparable across
    /// species — this is the quantity cross-species composition is built from.
    pub depth: f64,
    /// Total single-copy tags this cluster carries. `depth * n_markers` estimates the cluster's
    /// share of the sample's tag observations, which is what lets a caller normalize against
    /// **all** sequencing output rather than only the part that was resolved.
    pub n_markers: usize,
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
        let min_count = if p.adaptive_singleton { min_count_for(st.depth) } else { 2 };
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
    // Layer-1 selects the candidate clusters; Layer-2 quantifies them. The two are independent
    // so each can be A/B'd against the flat path on its own.
    let mut calls: Vec<StrainCall> = match p.layer1 {
        Layer1::Unique => profile_unique(db, counts, p),
        Layer1::Cst => match &db.tree {
            Some(tree) => descend_tree(tree, counts, p)
                .into_iter()
                .filter(|c| c.leaf < db.n_strains())
                .map(|c| StrainCall {
                    strain_index: c.leaf,
                    name: db.strain_names[c.leaf].clone(),
                    support: c.detected as f64,
                    coverage: c.coverage,
                    depth: c.depth,
                    n_markers: db.strain_markers[c.leaf].len(),
                    rel_abundance: 0.0,
                })
                .collect(),
            // A database built before the tree existed, or by `build` rather than `cluster`.
            // Degrade to the flat path rather than silently returning nothing.
            None => profile_unique(db, counts, p),
        },
    };

    if p.layer2 == Layer2::Enet && calls.len() > 1 {
        let idx: Vec<usize> = calls.iter().map(|c| c.strain_index).collect();
        let design = build_l2_design(db, &idx, counts);
        let selected = pre_scan(&design, 15, p.min_support_markers);
        if !selected.is_empty() {
            let w = l2_abundance(&design, &selected, p.enet_alpha, p.enet_l1_ratio);
            // Clusters the pre-scan dropped explained nothing beyond what the winners already
            // account for; keep only the selected ones, at their fitted depths.
            let keep: crate::fxhash::FxHashMap<usize, f64> = selected
                .iter()
                .zip(w.iter())
                .map(|(&col, &depth)| (design.clusters[col], depth))
                .collect();
            calls.retain(|c| keep.contains_key(&c.strain_index));
            for c in calls.iter_mut() {
                c.depth = keep[&c.strain_index];
            }
        }
    }

    normalize_by_depth(&mut calls);
    filter_by_abundance(&mut calls, p.min_rel_abundance);
    calls
}

/// The flat Layer-1: score each cluster independently on its own unique markers.
fn profile_unique(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<StrainCall> {
    let mut calls: Vec<StrainCall> = Vec::new();
    for j in 0..db.n_strains() {
        let st = panel_stats(db, counts, j);
        if st.panel == 0 {
            continue;
        }
        let min_count = if p.adaptive_singleton { min_count_for(st.depth) } else { 2 };
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
        // Depth–breadth consistency: reject a cluster whose markers are too DEEP for how FEW of
        // them were seen.
        //
        // This is the test that removes "shadow" clusters — the dominant false positive on real
        // panels. When the strain in the sample is not exactly any reference but sits between two
        // clusters, it carries all of cluster A's distinguishing loci and a fraction `f` of
        // cluster B's. Cluster B is then called on real reads at the *sample strain's* full
        // depth, but on only `f` of its panel. No coverage floor can catch this: measured on a
        // shadow scenario, B showed coverage 0.350 while a genuinely present B at 0.4x showed
        // 0.392 — indistinguishable. What differs is depth, 7.68x versus 0.44x.
        //
        // Under Poisson sampling a genuinely present cluster at depth λ must show breadth
        // `1 − e^(−λ)`. So `coverage / (1 − e^(−depth))` is ~1 for a real cluster at any depth,
        // and ~`f` for a shadow (its depth is `f × D`, large enough that the expected breadth is
        // ~1, while the observed breadth is only `f`). Swept on synthetic data: real clusters
        // scored 0.949–1.018 across depths 0.3x–20x, shadows scored 0.200/0.300/0.495/0.691/0.897
        // at f = 0.2/0.3/0.5/0.7/0.9. Real shadows carry small `f` — the S. epidermidis shadow in
        // MSA-1005 sat at ~7% of its true strain's abundance — so the default leaves wide margin
        // on both sides.
        //
        // Note the graceful degradation: as `f` → 1 a shadow becomes indistinguishable from a
        // genuine call, which is correct, because a strain carrying all of B's distinguishing
        // loci *is* evidence for B. The test is also inert below ~0.5x depth, where
        // `coverage ≈ 1 − e^(−depth)` holds for any cluster; a rare strain is never penalized.
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

    /// Regression: the abundance floor must not delete a correctly-estimated rare cluster.
    ///
    /// This is the interaction that collapsed recall on staggered mocks and high-host samples.
    /// StrainScan's 0.02 floor was calibrated against a *biased* depth estimator (median over
    /// detected markers, which pins a rare cluster near 1 read/tag regardless of how rare it is).
    /// Once the bias is removed the same floor cuts far deeper: here the rare cluster's true
    /// share is 1/(30+1) = 3.2%, and at a 30:0.5 depth ratio the correct answer is ~1.6% — which
    /// 0.02 would discard. Fixing an estimator without recalibrating the thresholds tuned to its
    /// bias is the failure mode this test exists to catch.
    #[test]
    fn abundance_floor_does_not_delete_correctly_estimated_rare_clusters() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![("dominant".into(), a.clone()), ("rare".into(), b.clone())]);

        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 30); // 30x
        }
        for &m in b.iter().take(400) {
            counts.insert(m, 1); // ~0.4x, 40% breadth
        }

        // Default params must keep both.
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

        // The old floor is still available, and still removes it — opt-in, not the default.
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

    /// A shadow cluster and a genuinely present rare cluster have the **same breadth** and differ
    /// only in depth, so only the depth–breadth consistency test can separate them.
    ///
    /// Shadow: the sample strain carries 30% of cluster B's distinguishing loci, so those markers
    /// appear at the sample strain's full 20x while the other 70% are absent — depth 6.0, breadth
    /// 0.30, and an expected breadth at that depth of ~1.0.
    /// Genuine: cluster B present at 0.33x — depth 0.33, breadth 0.33, expected breadth 0.28.
    #[test]
    fn shadow_clusters_are_rejected_but_genuine_rare_ones_are_kept() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B".into(), b.clone())]);

        let names = |p: &Params, counts: &MarkerCounts| -> Vec<String> {
            profile(&db, counts, p).iter().map(|c| c.name.clone()).collect()
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

        // --- shadow: 30% of B's panel at the true strain's depth
        let mut shadow = MarkerCounts::default();
        for &m in &a {
            shadow.insert(m, 20);
        }
        for &m in b.iter().take(300) {
            shadow.insert(m, 20);
        }
        assert_eq!(names(&loose, &shadow), vec!["A", "B"], "filter off: both called");
        assert_eq!(
            names(&default, &shadow),
            vec!["A"],
            "default must reject the shadow"
        );
        // Removing the shadow also repairs the true strain's abundance, which the shadow was
        // taking a 28% share of.
        let calls = profile(&db, &shadow, &default);
        assert!((calls[0].rel_abundance - 1.0).abs() < 1e-9);

        // --- genuine: B present at ~0.33x, i.e. the SAME breadth as the shadow
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

    /// **The reason the tree exists.** A leaf with too few markers of its own is invisible to the
    /// flat unique-marker algorithm, but the tree can accept it by pooling the markers of every
    /// ancestor whose sibling branch was never entered.
    ///
    /// Layout: 4 genomes forming ((A,B),(C,D)). A carries only 5 markers no one else has — below
    /// `min_support_markers` — but the A/B ancestor has 100 group-specific markers, and the root
    /// has 200 species-core markers. The sample contains A alone.
    ///
    /// Flat algorithm: A's own panel is 5 < 10, so A is never called.
    /// Tree: B and C/D are ruled out on their own sets, so nothing was entered on either sibling
    /// branch and A pools 5 + 100 + 200 = 305 markers, all observed.
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
        let ga = mk(&ab, 1000..1005); // only 5 exclusive markers
        let gb = mk(&ab, 1100..1200);
        let gc = mk(&cd, 1200..1300);
        let gd = mk(&cd, 1300..1400);
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> = [("A", &ga), ("B", &gb), ("C", &gc), ("D", &gd)]
            .into_iter()
            .map(|(n, g)| (n.to_string(), g.clone(), g.clone()))
            .collect();
        let cst = SpeciesCst::build(genomes, crate::cst::DEFAULT_SIMILARITY, false);
        assert_eq!(cst.n_clusters(), 4, "each genome should be its own cluster");

        // Sample: strain A only, at 20x.
        let mut counts = MarkerCounts::default();
        for &m in &ga {
            counts.insert(m, 20);
        }

        // --- flat algorithm: A has 5 unique markers, below the support floor -> missed
        let db = cst.cluster_db();
        let p = Params { min_rel_abundance: 0.0, ..Params::default() };
        let flat: Vec<String> = profile(&db, &counts, &p).iter().map(|c| c.name.clone()).collect();
        assert!(flat.is_empty(), "flat algorithm should miss the sparse leaf, got {flat:?}");

        // --- tree: pools ancestors whose siblings were never entered
        let tree = cst.build_tree();
        let calls = descend_tree(&tree, &counts, &p);
        assert_eq!(calls.len(), 1, "tree should call exactly the present leaf: {calls:?}");
        let call = &calls[0];
        assert_eq!(
            tree.leaves[call.leaf], vec![0],
            "the called leaf must be A (genome 0)"
        );
        assert_eq!(call.panel, 305, "pooled 5 own + 100 A/B-group + 200 root-core");
        assert!(call.path.len() == 3, "pooled leaf + 2 ancestors, got {:?}", call.path);
        assert!((call.coverage - 1.0).abs() < 1e-9);
        assert!((call.depth - 20.0).abs() < 0.5, "depth {}", call.depth);
    }

    /// Pooling must STOP once a sibling branch is entered: at that point the ancestor's markers
    /// are shared between both branches and attributing them to one would double-count.
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
        let genomes: Vec<(String, Vec<Marker>, Vec<Marker>)> = [("A", &ga), ("B", &gb), ("C", &gc), ("D", &gd)]
            .into_iter()
            .map(|(n, g)| (n.to_string(), g.clone(), g.clone()))
            .collect();
        let cst = SpeciesCst::build(genomes, crate::cst::DEFAULT_SIMILARITY, false);
        let tree = cst.build_tree();

        // BOTH A and B present -> the A/B ancestor is entered from both sides.
        let mut counts = MarkerCounts::default();
        for &m in ga.iter().chain(gb.iter()) {
            counts.insert(m, 20);
        }
        let p = Params { min_rel_abundance: 0.0, ..Params::default() };
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

    /// **The reason Layer-2 exists.** A cluster that is a strict subset of another has NO unique
    /// markers, so the flat algorithm cannot see it at all — `panel == 0` and it is skipped.
    /// The joint fit over shared markers recovers it exactly: the markers A and B share observe
    /// `w_A + w_B`, the ones only A carries observe `w_A`, and the two together give `w_B`.
    ///
    /// Sample: A at 10x and B at 5x, B's marker set a strict subset of A's.
    #[test]
    fn joint_fit_recovers_a_subset_cluster_the_flat_path_cannot_see() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let b: Vec<Marker> = (10_000..10_500).collect(); // strict subset of A
        let db = StrainDb::build(vec![("A".into(), a.clone()), ("B_subset".into(), b.clone())]);
        assert_eq!(db.unique_marker_count(1), 0, "B has no unique markers by construction");

        let mut counts = MarkerCounts::default();
        for &m in &a {
            // shared rows see both strains; A-only rows see A alone
            counts.insert(m, if b.contains(&m) { 15 } else { 10 });
        }

        // Flat path: B is invisible.
        let p = Params { min_rel_abundance: 0.0, ..Params::default() };
        let flat: Vec<String> = profile(&db, &counts, &p).iter().map(|c| c.name.clone()).collect();
        assert_eq!(flat, vec!["A".to_string()], "flat path cannot see the subset cluster");

        // Joint fit: both recovered, at the right depths.
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

    /// The pre-scan must consume the winner's markers so a near-duplicate cluster is not selected
    /// on the same evidence. With unique-only markers this loop is a no-op — the sets are
    /// disjoint, so consuming one leaves the others untouched — which is why it only becomes
    /// meaningful once shared markers are in the matrix.
    #[test]
    fn pre_scan_consumes_the_winners_markers() {
        let a: Vec<Marker> = (10_000..11_000).collect();
        let dup: Vec<Marker> = (10_000..10_990).collect(); // 99% the same as A
        let far: Vec<Marker> = (20_000..21_000).collect();
        let db = StrainDb::build(vec![
            ("A".into(), a.clone()),
            ("A_dup".into(), dup),
            ("Far".into(), far.clone()),
        ]);
        let mut counts = MarkerCounts::default();
        for &m in &a {
            counts.insert(m, 20); // only A is present
        }

        let design = build_l2_design(&db, &[0, 1, 2], &counts);
        let chosen = pre_scan(&design, 15, 10);
        assert_eq!(chosen.first(), Some(&0), "A explains the most, picked first");
        assert!(
            !chosen.contains(&1),
            "the 99%-duplicate has almost nothing left to explain once A's markers are consumed"
        );
        assert!(!chosen.contains(&2), "the absent cluster explains nothing");
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

// ===== Layer-1: Cluster Search Tree descent (StrainScan port) ==============

use crate::cst::Cst;

/// Minimum markers for a node's own set to be worth testing. StrainScan uses 1000 k-mers, but it
/// stores both strands as separate ids (so 500 canonical) and works on a ~76-116x denser marker
/// space; the scale-equivalent here is single digits. 10 is that floor rounded up to the value
/// already used for cluster support, so the two gates stay on one scale.
pub const MIN_NODE_MARKERS: usize = 10;

/// One leaf accepted by the tree descent.
#[derive(Debug, Clone)]
pub struct TreeCall {
    pub leaf: usize,
    /// Markers pooled along the unique path (the leaf's own plus attributable ancestors').
    pub panel: usize,
    pub detected: usize,
    pub coverage: f64,
    pub depth: f64,
    /// Nodes whose markers were pooled — leaf first, then ancestors.
    pub path: Vec<usize>,
}

/// Evidence for one marker set.
fn set_evidence(markers: &[Marker], counts: &MarkerCounts) -> (usize, usize, f64) {
    let panel = markers.len();
    if panel == 0 {
        return (0, 0, 0.0);
    }
    let mut obs: Vec<u32> = markers
        .iter()
        .map(|m| counts.get(m).copied().unwrap_or(0))
        .collect();
    obs.sort_unstable();
    let detected = panel - obs.partition_point(|&c| c < 1);
    let depth = if detected == 0 {
        0.0
    } else {
        let cap_idx = panel.saturating_sub(1 + detected / TRIM_FRACTION);
        let cap = obs[cap_idx] as u64;
        obs.iter().map(|&c| (c as u64).min(cap)).sum::<u64>() as f64 / panel as f64
    };
    (panel, detected, depth)
}

/// Descend the Cluster Search Tree, returning the leaves it accepts.
///
/// Two mechanisms matter, and they are the entire reason to build a tree:
///
/// **Pruning.** A node whose own marker set is informative but unobserved rules out its whole
/// subtree in one test, instead of scoring every leaf independently.
///
/// **Unique-path pooling.** A leaf is accepted on the union of its own markers *and* those of
/// every ancestor whose sibling branch was never entered. If the descent went left at a node and
/// never right, that node's group-specific markers are attributable to the left subtree, so a
/// leaf with too few markers of its own can borrow them. This is what the flat unique-marker
/// algorithm has no way to do: it sees only the leaf's own set and calls the leaf undetectable.
/// Once a sibling *is* entered, the ancestor's markers become ambiguous between the two branches
/// and pooling stops there.
pub fn descend_tree(cst: &Cst, counts: &MarkerCounts, p: &Params) -> Vec<TreeCall> {
    if cst.n_leaves() == 0 {
        return Vec::new();
    }
    if cst.n_leaves() == 1 {
        // Degenerate tree: the single leaf is the root; test it directly.
        let ms: Vec<Marker> = cst.node_markers[0].iter().copied().collect();
        let (panel, detected, depth) = set_evidence(&ms, counts);
        if panel > 0 && detected >= p.min_support_markers {
            let coverage = detected as f64 / panel as f64;
            if coverage >= p.min_coverage {
                return vec![TreeCall { leaf: 0, panel, detected, coverage, depth, path: vec![0] }];
            }
        }
        return Vec::new();
    }

    // A node "fires" when its own marker set is informative AND observed.
    let fires = |v: usize| -> bool {
        let ms: Vec<Marker> = cst.node_markers[v].iter().copied().collect();
        if ms.len() < MIN_NODE_MARKERS {
            return false; // uninformative: cannot rule the subtree in or out
        }
        let (panel, detected, depth) = set_evidence(&ms, counts);
        let min_count = if p.adaptive_singleton { min_count_for(depth) } else { 2 };
        let support = if min_count >= 2 {
            let mut n = 0;
            for m in &ms {
                if counts.get(m).copied().unwrap_or(0) >= 2 {
                    n += 1;
                }
            }
            n
        } else {
            detected
        };
        support >= p.min_support_markers && (detected as f64 / panel as f64) >= p.min_coverage
    };
    let informative = |v: usize| cst.node_markers[v].len() >= MIN_NODE_MARKERS;

    // Descend, recording which nodes were entered.
    let mut entered: Vec<bool> = vec![false; cst.n_nodes()];
    let mut reached: Vec<usize> = Vec::new();
    let mut stack = vec![cst.root];
    entered[cst.root] = true;
    while let Some(v) = stack.pop() {
        if cst.is_leaf(v) {
            reached.push(v);
            continue;
        }
        let (a, b) = cst.children[v].expect("internal node has children");
        for c in [a, b] {
            // Enter a child if it has its own evidence, or if it cannot be tested at all.
            // An uninformative child is not evidence of absence, so we descend and let a
            // deeper node decide.
            if !informative(c) || fires(c) {
                entered[c] = true;
                stack.push(c);
            }
        }
    }

    // Accept reached leaves on pooled evidence along the unique path.
    let mut out = Vec::new();
    for &leaf in &reached {
        let mut path = vec![leaf];
        let mut pooled: Vec<Marker> = cst.node_markers[leaf].iter().copied().collect();
        let mut v = leaf;
        while let Some(par) = cst.parent[v] {
            match cst.sibling(v) {
                // The sibling branch was never entered, so this ancestor's group-specific
                // markers belong to our side and can be pooled.
                Some(s) if !entered[s] => {
                    pooled.extend(cst.node_markers[par].iter().copied());
                    path.push(par);
                    v = par;
                }
                // Sibling entered: the ancestor's markers are shared between both branches and
                // attributing them here would double-count. Stop.
                _ => break,
            }
        }
        let (panel, detected, depth) = set_evidence(&pooled, counts);
        if panel == 0 {
            continue;
        }
        let min_count = if p.adaptive_singleton { min_count_for(depth) } else { 2 };
        let support = if min_count >= 2 {
            pooled
                .iter()
                .filter(|m| counts.get(m).copied().unwrap_or(0) >= 2)
                .count()
        } else {
            detected
        };
        if support < p.min_support_markers {
            continue;
        }
        let coverage = detected as f64 / panel as f64;
        if coverage < p.min_coverage {
            continue;
        }
        let expected = detectable_fraction(depth);
        if expected > 0.0 && coverage / expected < p.min_consistency {
            continue;
        }
        out.push(TreeCall { leaf, panel, detected, coverage, depth, path });
    }
    out
}

// ===== Layer-2: shared-marker deconvolution (StrainScan port) ==============
//
// The flat algorithm scores each cluster on markers **no other cluster carries**, which on a
// low-diversity species is a small minority of the panel. StrainScan instead fits all co-present
// clusters jointly over the *shared* markers too: a marker carried by clusters {A,C} constrains
// `w_A + w_C` against its observed count, which is information the unique-only path throws away.
//
// Three pieces, in StrainScan's order:
//   1. `L2Design`      — the marker × cluster incidence matrix, shared markers included.
//   2. `pre_scan`      — greedy selection by *residually* covered markers, consuming each
//                        winner's markers so the next candidate is judged on what is left.
//   3. `nonneg_elastic_net` — joint abundance over the selected columns.

/// Marker × cluster incidence for one species' co-detected clusters.
#[derive(Debug, Clone)]
pub struct L2Design {
    /// Row order: the markers used, each carried by 1..n of the candidate clusters.
    pub markers: Vec<Marker>,
    /// Column order: indices into the `StrainDb`.
    pub clusters: Vec<usize>,
    /// `cols[j][i] == 1.0` iff cluster `clusters[j]` carries `markers[i]`.
    pub cols: Vec<Vec<f64>>,
    /// Observed count per row.
    pub y: Vec<f64>,
}

impl L2Design {
    pub fn n_rows(&self) -> usize {
        self.markers.len()
    }
    /// Fraction of rows carried by more than one candidate — the data the unique-only path
    /// discards. If this is ~0 the clusters share nothing and a joint fit cannot beat
    /// per-cluster means.
    pub fn shared_fraction(&self) -> f64 {
        if self.markers.is_empty() {
            return 0.0;
        }
        let shared = (0..self.markers.len())
            .filter(|&i| self.cols.iter().filter(|c| c[i] > 0.0).count() > 1)
            .count();
        shared as f64 / self.markers.len() as f64
    }
}

/// Build the design matrix over `candidates` within one species database.
///
/// A marker enters as a row when at least one candidate carries it — **including** markers all
/// of them carry.
///
/// This is a deliberate deviation from StrainScan, which drops the all-carried rows. It is right
/// there and wrong here, because the two setups differ in what is already known. StrainScan's
/// Layer-2 candidates are strains *within one cluster* whose total depth Layer-1 has already
/// pinned down, so a row carried by every candidate adds no information about the split. Our
/// candidates are co-detected *clusters* with no separately determined total, so those rows are
/// the only thing constraining the sum: with A ⊃ B, the markers they share observe `w_A + w_B`
/// and the ones only A carries observe `w_A`, and it takes both to recover `w_B`. Dropping them
/// would break exactly the case this matrix exists for.
pub fn build_l2_design(db: &StrainDb, candidates: &[usize], counts: &MarkerCounts) -> L2Design {
    let mut markers: Vec<Marker> = Vec::new();
    let mut seen: crate::fxhash::FxHashSet<Marker> = crate::fxhash::FxHashSet::default();
    for &j in candidates {
        for &m in &db.strain_markers[j] {
            if seen.insert(m) {
                markers.push(m);
            }
        }
    }
    markers.sort_unstable();

    let cols: Vec<Vec<f64>> = candidates
        .iter()
        .map(|&j| {
            markers
                .iter()
                .map(|m| if db.strain_markers[j].contains(m) { 1.0 } else { 0.0 })
                .collect()
        })
        .collect();
    let y: Vec<f64> = markers
        .iter()
        .map(|m| counts.get(m).copied().unwrap_or(0) as f64)
        .collect();
    L2Design { markers, clusters: candidates.to_vec(), cols, y }
}

/// StrainScan's iterative pre-scan: greedily pick the cluster explaining the most **not yet
/// explained** markers, then consume its markers so later candidates are judged on the residual.
///
/// This is the step the unique-only path cannot have: cluster-unique sets are disjoint by
/// construction, so consuming one leaves the others untouched and the iteration is a no-op. It
/// only bites once shared markers are in the matrix.
///
/// Returns column indices into `design.clusters`, in selection order.
pub fn pre_scan(design: &L2Design, max_iter: usize, min_new_markers: usize) -> Vec<usize> {
    let n_rows = design.n_rows();
    let n_cols = design.cols.len();
    if n_rows == 0 || n_cols == 0 {
        return Vec::new();
    }
    let mut used = vec![false; n_rows];
    let mut chosen: Vec<usize> = Vec::new();
    let mut taken = vec![false; n_cols];

    for _ in 0..max_iter.min(n_cols) {
        let mut best = (0usize, 0usize); // (residual support, column)
        for (j, is_taken) in taken.iter().enumerate() {
            if *is_taken {
                continue;
            }
            // Score by how many *residual* markers this cluster explains at count >= 1.
            let score = (0..n_rows)
                .filter(|&i| !used[i] && design.cols[j][i] > 0.0 && design.y[i] >= 1.0)
                .count();
            if score > best.0 {
                best = (score, j);
            }
        }
        if best.0 < min_new_markers {
            break;
        }
        let j = best.1;
        chosen.push(j);
        taken[j] = true;
        // Consume this cluster's markers unconditionally — this is what makes the loop
        // terminate and what stops a near-duplicate cluster being selected on the same
        // evidence as the winner.
        for (i, u) in used.iter_mut().enumerate() {
            if design.cols[j][i] > 0.0 {
                *u = true;
            }
        }
    }
    chosen
}

/// Joint abundance for the selected clusters, via the non-negative Elastic Net.
///
/// Returns one depth per entry of `selected` (indices into `design.clusters`), in the same order.
pub fn l2_abundance(
    design: &L2Design,
    selected: &[usize],
    alpha: f64,
    l1_ratio: f64,
) -> Vec<f64> {
    if selected.is_empty() || design.n_rows() == 0 {
        return vec![0.0; selected.len()];
    }
    let cols: Vec<Vec<f64>> = selected.iter().map(|&j| design.cols[j].clone()).collect();
    nonneg_elastic_net(&cols, &design.y, alpha, l1_ratio, 2000, 1e-8)
}
