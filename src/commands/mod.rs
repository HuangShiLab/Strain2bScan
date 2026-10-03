//! Subcommand implementations for the strain2bscan CLI.

use strain2bscan::identify::{Layer1, Layer2, Params};

pub mod batch;
pub mod build;
pub mod cluster;
pub mod cst_demo;
pub mod demo;
pub mod diagnose_tree;
pub mod evaluate;
pub mod info;
pub mod multi_profile;
pub mod profile;

#[allow(clippy::too_many_arguments)]
pub(crate) fn parse_params(
    min_support: Option<usize>,
    min_coverage: Option<f64>,
    min_abundance: Option<f64>,
    trace_gap: Option<f64>,
    trace_floor: Option<f64>,
    layer1: Option<&str>,
    layer2: Option<&str>,
    enet_alpha: Option<f64>,
    min_consistency: Option<f64>,
    fixed_gate: bool,
    no_adaptive_singleton: bool,
    no_adaptive_floor: bool,
) -> Result<Params, String> {
    let mut p = Params::default();
    if let Some(v) = min_support {
        p.min_support_markers = v;
    }
    if let Some(v) = min_coverage {
        p.min_coverage = v;
    }
    if let Some(v) = min_abundance {
        p.min_rel_abundance = v;
    }
    if let Some(v) = trace_gap {
        p.trace_gap = v;
    }
    if let Some(v) = trace_floor {
        p.trace_floor = v;
    }
    match layer1 {
        None | Some("auto") => {}
        Some("unique") => p.layer1 = Layer1::Unique,
        Some("cst") => p.layer1 = Layer1::Cst,
        Some(x) => return Err(format!("bad --layer1 {x} (want auto|unique|cst)")),
    }
    match layer2 {
        None | Some("depth") => {}
        Some("enet") => p.layer2 = Layer2::Enet,
        Some(x) => return Err(format!("bad --layer2 {x} (want depth|enet)")),
    }
    if let Some(v) = enet_alpha {
        p.enet_alpha = v;
    }
    if let Some(v) = min_consistency {
        p.min_consistency = v;
    }
    if fixed_gate {
        p.adaptive_singleton = false;
        p.adaptive_floor = false;
    }
    if no_adaptive_singleton {
        p.adaptive_singleton = false;
    }
    if no_adaptive_floor {
        p.adaptive_floor = false;
    }
    Ok(p)
}
