//! Presence detection by unique markers.

use crate::db::StrainDb;
use crate::depth::support_count;
use crate::identify::Params;
use crate::markers::MarkerCounts;

/// Detect present clusters/strains by their **unique** markers only.
///
/// Returns `(cluster_index, supporting_marker_count)`.
pub fn detect_present(db: &StrainDb, counts: &MarkerCounts, p: &Params) -> Vec<(usize, f64)> {
    let mut out = Vec::new();
    for j in 0..db.n_strains() {
        let st = panel_stats(db, counts, j);
        let support = support_count(&st, p.adaptive_singleton);
        if support >= p.min_support_markers {
            out.push((j, support as f64));
        }
    }
    out
}

// Re-export the panel-stat helper used here so callers can stay in `identify::detect`.
pub use crate::depth::panel_stats;
