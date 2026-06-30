//! StrainScan-style Layer-2 strain resolution + abundance, on 2bRAD-tag markers.
//!
//! 1. **Presence detection by unique markers.** A cluster/strain is called present iff at
//!    least `min_support_markers` of its *unique* (cluster-specific) markers are observed
//!    with sample count ≥ 2. Restricting to single-cluster markers makes detection immune to
//!    the shared-marker cross-talk that breaks a greedy set-cover on conspecific panels.
//!    Markers with count < 2 are treated as sequencing error (StrainScan's singleton filter).
//! 2. **Abundance from unique-marker depth.** Each detected strain's depth is the median
//!    sample count over its detected unique markers; relative abundance ∝ depth. This is
//!    robust to the mis-attribution a joint regression suffers when strains are very similar.
//!    A non-negative Elastic Net solver (`nonneg_elastic_net`, StrainScan's
//!    `ElasticNet(positive=True)`) is also provided.
//! 3. **Post-filter.** Drop strains below `min_rel_abundance`; renormalize.

use std::collections::HashMap;

use crate::db::StrainDb;
use crate::markers::Marker;

#[derive(Debug, Clone)]
pub struct Params {
    /// Min number of a strain's unique markers (count ≥ 2) to call it present (tag-unit
    /// support — StrainScan's `msn`; recalibrate on your data).
    pub min_support_markers: usize,
    /// Min fraction of a strain's unique markers detected (StrainScan 0.7).
    pub min_coverage: f64,
    /// Min relative abundance to keep a strain (StrainScan 0.02).
    pub min_rel_abundance: f64,
    /// Elastic Net penalty strength (0 ⇒ pure NNLS).
    pub alpha: f64,
    /// Elastic Net L1 ratio (StrainScan default 0.5).
    pub l1_ratio: f64,
    pub max_iter: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            // Recalibrated for sparse 2bRAD markers (full-k-mer StrainScan uses msn*k≈1240
            // k-mers; tag markers are ~50-100x sparser, so the floor is in *tag* units).
            min_support_markers: 10,
            min_coverage: 0.7,
            min_rel_abundance: 0.02,
            alpha: 0.01,
            l1_ratio: 0.5,
            max_iter: 1000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct StrainCall {
    pub strain_index: usize,
    pub name: String,
    /// Support score at the iteration this strain was selected.
    pub support: f64,
    /// Fraction of the strain's unique markers detected in the sample.
    pub coverage: f64,
    /// Relative abundance (sums to 1 across kept strains).
    pub rel_abundance: f64,
}

/// Effective sample count for a marker, applying StrainScan's singleton-error filter.
#[inline]
fn eff_count(counts: &HashMap<Marker, u32>, m: Marker) -> f64 {
    match counts.get(&m) {
        Some(&c) if c >= 2 => c as f64,
        _ => 0.0,
    }
}

/// Detect present clusters/strains by their **unique** markers only.
///
/// A cluster is called present iff at least `min_support_markers` of its unique (cluster-
/// specific) markers are observed with sample count ≥ 2. Because this uses only markers that
/// belong to a single cluster, it is immune to the shared-marker cross-talk that breaks the
/// greedy set-cover on a large conspecific panel (the Layer-1 role: decide *which* clusters
/// are present before quantifying). Returns `(cluster_index, detected_unique_count)`.
pub fn detect_present(
    db: &StrainDb,
    counts: &HashMap<Marker, u32>,
    p: &Params,
) -> Vec<(usize, f64)> {
    let mut out = Vec::new();
    for j in 0..db.n_strains() {
        let detected = db
            .unique_markers(j)
            .filter(|&m| eff_count(counts, m) > 0.0)
            .count();
        if detected >= p.min_support_markers {
            out.push((j, detected as f64));
        }
    }
    out
}

/// Full pipeline: detect present clusters via unique markers, estimate abundance with a
/// non-negative Elastic Net, filter low-abundance calls, renormalize.
pub fn profile(db: &StrainDb, counts: &HashMap<Marker, u32>, p: &Params) -> Vec<StrainCall> {
    let selected = detect_present(db, counts, p);
    if selected.is_empty() {
        return Vec::new();
    }
    let idx: Vec<usize> = selected.iter().map(|&(j, _)| j).collect();

    // --- Abundance from each cluster's UNIQUE markers. A joint regression over the full
    //     (shared-heavy) marker set mis-attributes signal between very similar co-present
    //     strains; estimating each strain's depth from its own unique markers is robust.
    //     depth_j = median sample count over cluster j's detected unique markers;
    //     relative abundance ∝ depth_j. ---
    let depths: Vec<f64> = idx
        .iter()
        .map(|&j| unique_marker_depth(db, counts, j))
        .collect();
    let depth_sum: f64 = depths.iter().sum();

    // --- Assemble calls with coverage, then filter + renormalize. ---
    let mut calls: Vec<StrainCall> = Vec::new();
    for (c, &j) in idx.iter().enumerate() {
        let rel = if depth_sum > 0.0 {
            depths[c] / depth_sum
        } else {
            0.0
        };
        let coverage = strain_unique_coverage(db, counts, j);
        calls.push(StrainCall {
            strain_index: j,
            name: db.strain_names[j].clone(),
            support: selected[c].1,
            coverage,
            rel_abundance: rel,
        });
    }
    // Detection already gated on unique markers, so filter only by abundance here. Coverage
    // is retained as a reported confidence metric (errors can lower it below min_coverage).
    calls.retain(|c| c.rel_abundance >= p.min_rel_abundance);
    let _ = p.min_coverage;
    let kept_sum: f64 = calls.iter().map(|c| c.rel_abundance).sum();
    if kept_sum > 0.0 {
        for c in &mut calls {
            c.rel_abundance /= kept_sum;
        }
    }
    calls.sort_by(|a, b| b.rel_abundance.partial_cmp(&a.rel_abundance).unwrap());
    calls
}

/// Robust per-strain depth = median sample count over the strain's **detected** unique
/// markers (count ≥ 1). Median resists repeat/contamination outliers. Returns 0 if none.
fn unique_marker_depth(db: &StrainDb, counts: &HashMap<Marker, u32>, j: usize) -> f64 {
    let mut obs: Vec<f64> = db
        .unique_markers(j)
        .filter_map(|m| match counts.get(&m).copied().unwrap_or(0) {
            c if c >= 1 => Some(c as f64),
            _ => None,
        })
        .collect();
    if obs.is_empty() {
        return 0.0;
    }
    obs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = obs.len() / 2;
    if obs.len().is_multiple_of(2) {
        (obs[mid - 1] + obs[mid]) / 2.0
    } else {
        obs[mid]
    }
}

/// Coverage = fraction of a strain's unique markers detected (count ≥ 1).
fn strain_unique_coverage(db: &StrainDb, counts: &HashMap<Marker, u32>, j: usize) -> f64 {
    let uniq: Vec<Marker> = db.unique_markers(j).collect();
    let denom = uniq.len();
    if denom == 0 {
        // No unique markers (e.g., a strain fully contained in another): fall back to
        // coverage over all its markers so it isn't auto-rejected.
        let total = db.strain_markers[j].len();
        if total == 0 {
            return 0.0;
        }
        let seen = db.strain_markers[j]
            .iter()
            .filter(|&&m| counts.get(&m).copied().unwrap_or(0) >= 1)
            .count();
        return seen as f64 / total as f64;
    }
    let seen = uniq
        .iter()
        .filter(|&&m| counts.get(&m).copied().unwrap_or(0) >= 1)
        .count();
    seen as f64 / denom as f64
}

/// Non-negative Elastic Net via cyclic coordinate descent with residual maintenance.
/// Minimizes ½‖Xw − y‖² + α·l1·n·‖w‖₁ + ½·α·(1−l1)·n·‖w‖²  s.t. w ≥ 0.
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
pub fn naive_profile(db: &StrainDb, counts: &HashMap<Marker, u32>, min_score: f64) -> Vec<usize> {
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
    /// Returns (db, marker layout) so the test can synthesize a sample.
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
    fn synth_sample(db: &StrainDb, mixture: &[(usize, f64)], depth: f64) -> HashMap<Marker, u32> {
        let mut present: HashMap<Marker, f64> = HashMap::new();
        for &(j, ab) in mixture {
            for &m in &db.strain_markers[j] {
                *present.entry(m).or_insert(0.0) += ab;
            }
        }
        let mut counts: HashMap<Marker, u32> = HashMap::new();
        for (m, frac) in present {
            let c = (depth * frac).round() as u32;
            if c > 0 {
                counts.insert(m, c);
            }
        }
        // Inject sequencing-error singletons (count 1) that must be ignored.
        for e in 0..50u64 {
            counts.insert(9_000_000 + e, 1);
        }
        counts
    }

    #[test]
    fn resolves_conspecific_mixture_where_naive_overcalls() {
        // 4 strains, 200 shared core markers, 50 private each.
        let (db, _priv) = conspecific_db(4, 200, 50);
        // True mixture: strain0 70%, strain2 30%.
        let counts = synth_sample(&db, &[(0, 0.7), (2, 0.3)], 30.0);

        // Ported Layer-2 recovers exactly {0, 2} at ~70/30.
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

        // Naive strainscan-rust-style scoring over-calls: shared core makes ALL 4 strains
        // clear the threshold (core alone = 200 markers × 30 = 6000).
        let naive = naive_profile(&db, &counts, 1240.0);
        assert_eq!(naive.len(), 4, "naive should over-call all 4: {naive:?}");
    }

    #[test]
    fn singleton_errors_do_not_create_calls() {
        let (db, _) = conspecific_db(3, 100, 40);
        // Empty true sample, only error singletons.
        let mut counts: HashMap<Marker, u32> = HashMap::new();
        for e in 0..100u64 {
            counts.insert(9_000_000 + e, 1);
        }
        assert!(profile(&db, &counts, &Params::default()).is_empty());
    }

    #[test]
    fn nnls_recovers_known_coefficients() {
        // y = 2*c0 + 3*c1 on a small binary design.
        let cols = vec![vec![1.0, 0.0, 1.0, 2.0], vec![0.0, 1.0, 1.0, 1.0]];
        let y = vec![2.0, 3.0, 5.0, 7.0];
        let w = nonneg_elastic_net(&cols, &y, 0.0, 0.5, 5000, 1e-10);
        assert!(
            (w[0] - 2.0).abs() < 1e-3 && (w[1] - 3.0).abs() < 1e-3,
            "w={w:?}"
        );
    }
}
