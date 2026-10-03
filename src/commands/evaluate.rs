//! Evaluate predictions against a truth table.

use crate::cli::EvaluateArgs;
use strain2bscan::bench::{evaluate, parse_abundance_tsv};

pub fn run(args: &EvaluateArgs) -> Result<(), String> {
    let pred_text = std::fs::read_to_string(&args.pred).map_err(|e| e.to_string())?;
    let truth_text = std::fs::read_to_string(&args.truth).map_err(|e| e.to_string())?;
    let pred = parse_abundance_tsv(&pred_text);
    let truth = parse_abundance_tsv(&truth_text);
    let m = evaluate(&pred, &truth, args.present);
    println!(
        "TP={} FP={} FN={}  precision={:.3} recall={:.3} F1={:.3}  L1={:.3} Bray-Curtis={:.3}",
        m.tp, m.fp, m.fn_, m.precision, m.recall, m.f1, m.l1, m.bray_curtis
    );
    Ok(())
}
