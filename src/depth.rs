//! Depth estimation and marker-evidence helpers for strain profiling.

use crate::db::StrainDb;
use crate::markers::{Marker, MarkerCounts};

/// Depth at or above which a genuine marker is essentially never observed exactly once, so
/// `count == 1` can safely be attributed to sequencing error (StrainScan's singleton rule).
pub const SINGLETON_SAFE_DEPTH: f64 = 3.0;

/// Reciprocal of the fraction of **non-zero** observations winsorized before averaging, to keep
/// collapsed repeats and contamination from inflating the depth estimate (top 1%).
const TRIM_FRACTION: usize = 100;

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
pub struct PanelStats {
    /// Panel size (unique markers, or all markers when the cluster has no unique ones).
    pub panel: usize,
    /// Markers with count >= 1.
    pub detected1: usize,
    /// Markers with count >= 2.
    pub detected2: usize,
    /// Zero-inclusive trimmed mean count — the absolute depth estimate.
    pub depth: f64,
}

/// Compute panel size, detected counts (>=1 and >=2), and zero-inclusive winsorized depth over an
/// arbitrary marker set.
///
/// This is the single implementation of the depth estimator: mean count over the whole panel
/// (zeros included), with the top 1% of *non-zero* observations winsorized down to the 99th
/// percentile. Keeping it in one place removes the divergence risk between the flat unique-marker
/// path ([`panel_stats`]) and the tree-descent path ([`set_evidence`](crate::tree::set_evidence)).
pub fn marker_panel_evidence(markers: &[Marker], counts: &MarkerCounts) -> PanelStats {
    let panel = markers.len();
    if panel == 0 {
        return PanelStats::default();
    }

    // One pass: count markers seen once / twice, and collect only the non-zero counts for the
    // winsorized depth estimate. Avoiding a full sort on the whole panel is a noticeable win when
    // `panel_stats` is called once per cluster.
    let mut detected1 = 0usize;
    let mut detected2 = 0usize;
    let mut nonzero: Vec<u32> = Vec::new();
    for &m in markers {
        let c = counts.get(&m).copied().unwrap_or(0);
        if c >= 1 {
            detected1 += 1;
            nonzero.push(c);
        }
        if c >= 2 {
            detected2 += 1;
        }
    }

    // Depth = mean count over the WHOLE panel (zeros included — they are the evidence that the
    // strain is rare), with the top 1% of *non-zero* observations winsorized down to the 99th
    // percentile so collapsed repeats and contamination cannot inflate it.
    let depth = if detected1 == 0 {
        0.0
    } else {
        let trim = detected1 / TRIM_FRACTION;
        let cap = if trim == 0 {
            *nonzero.iter().max().unwrap() as u64
        } else {
            let k = detected1 - trim - 1;
            *nonzero.select_nth_unstable(k).1 as u64
        };
        let sum: u64 = nonzero.iter().map(|&c| (c as u64).min(cap)).sum();
        sum as f64 / panel as f64
    };

    PanelStats {
        panel,
        detected1,
        detected2,
        depth,
    }
}

/// Compute [`PanelStats`] over cluster `j`'s **unique** markers.
pub fn panel_stats(db: &StrainDb, counts: &MarkerCounts, j: usize) -> PanelStats {
    marker_panel_evidence(db.unique_markers(j), counts)
}

/// Robust per-strain absolute depth: the **zero-inclusive** trimmed mean count over the
/// strain's unique-marker panel.
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

/// Number of markers that meet the count threshold.
///
/// This is the single place where the singleton policy is applied, so the tree descent, the flat
/// unique-marker path, and the clade fallback cannot drift out of agreement.
pub fn support_count(ev: &PanelStats, adaptive_singleton: bool) -> usize {
    let min_count = if adaptive_singleton {
        min_count_for(ev.depth)
    } else {
        2
    };
    if min_count >= 2 {
        ev.detected2
    } else {
        ev.detected1
    }
}
