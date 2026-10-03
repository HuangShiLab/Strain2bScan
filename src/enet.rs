//! Layer-2: shared-marker deconvolution via non-negative Elastic Net.

use crate::db::StrainDb;
use crate::depth::marker_panel_evidence;
use crate::identify::Params;
use crate::markers::{Marker, MarkerCounts};

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
    let mut r = y.to_vec();
    let col_sq: Vec<f64> = cols.iter().map(|c| c.iter().map(|v| v * v).sum()).collect();
    let l1 = alpha * l1_ratio * n as f64;
    let l2 = alpha * (1.0 - l1_ratio) * n as f64;

    for _ in 0..max_iter {
        let mut max_dw = 0.0_f64;
        for j in 0..k {
            if col_sq[j] == 0.0 {
                continue;
            }
            let mut rho = col_sq[j] * w[j];
            for i in 0..n {
                rho += cols[j][i] * r[i];
            }
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

/// Marker × cluster incidence for one species' co-detected clusters.
#[derive(Debug, Clone)]
pub struct L2Design {
    /// Row order: the markers used.
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

/// Minimum share of the total fitted depth for a cluster with **no unique evidence** to be
/// reported.
pub const MIN_SUBSET_SHARE: f64 = 0.02;

/// Clusters that Layer-1 **structurally cannot see**, but that the shared-marker fit can resolve.
pub fn subset_candidates(
    db: &StrainDb,
    counts: &MarkerCounts,
    called: &[usize],
    p: &Params,
) -> Vec<usize> {
    if called.is_empty() {
        return Vec::new();
    }
    let is_called = |j: usize| called.contains(&j);

    let mut out = Vec::new();
    for j in 0..db.n_strains() {
        if is_called(j) || db.unique_marker_count(j) >= p.min_support_markers {
            continue;
        }
        let ms = &db.strain_markers[j];
        if ms.is_empty() {
            continue;
        }
        let outside = called
            .iter()
            .map(|&k| {
                let other = &db.strain_markers[k];
                ms.iter().filter(|m| !other.contains(m)).count()
            })
            .min()
            .unwrap_or(usize::MAX);
        if outside >= p.min_support_markers {
            continue;
        }
        let v: Vec<Marker> = ms.iter().copied().collect();
        let ev = marker_panel_evidence(&v, counts);
        if ev.panel == 0 || (ev.detected1 as f64 / ev.panel as f64) < p.min_coverage {
            continue;
        }
        out.push(j);
    }
    out
}

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
                .map(|m| {
                    if db.strain_markers[j].contains(m) {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let y: Vec<f64> = markers
        .iter()
        .map(|m| counts.get(m).copied().unwrap_or(0) as f64)
        .collect();
    L2Design {
        markers,
        clusters: candidates.to_vec(),
        cols,
        y,
    }
}

/// StrainScan's iterative pre-scan: greedily pick the cluster explaining the most residual
/// markers, then consume its markers.
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
        let mut best = (0usize, 0usize);
        for (j, is_taken) in taken.iter().enumerate() {
            if *is_taken {
                continue;
            }
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
        for (i, u) in used.iter_mut().enumerate() {
            if design.cols[j][i] > 0.0 {
                *u = true;
            }
        }
    }
    chosen
}

/// Joint abundance for the selected clusters, via the non-negative Elastic Net.
pub fn l2_abundance(design: &L2Design, selected: &[usize], alpha: f64, l1_ratio: f64) -> Vec<f64> {
    if selected.is_empty() || design.n_rows() == 0 {
        return vec![0.0; selected.len()];
    }
    let cols: Vec<Vec<f64>> = selected.iter().map(|&j| design.cols[j].clone()).collect();
    nonneg_elastic_net(&cols, &design.y, alpha, l1_ratio, 2000, 1e-8)
}
