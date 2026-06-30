//! Benchmark metrics for strain-level profiling against MockMetagenomes4Benchmark.
//!
//! Given predicted relative abundances and a ground-truth composition (both as
//! `name -> abundance`), compute the standard accuracy metrics used to evaluate
//! strain/cluster profilers. Detection is judged at a presence threshold; abundance error
//! uses L1 and Bray–Curtis over the union of predicted and truth labels.

use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct Metrics {
    pub tp: usize,
    pub fp: usize,
    pub fn_: usize,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
    /// L1 distance between abundance vectors over the label union (0=identical, 2=disjoint).
    pub l1: f64,
    /// Bray–Curtis dissimilarity (0=identical, 1=disjoint).
    pub bray_curtis: f64,
}

fn presence(map: &HashMap<String, f64>, thresh: f64) -> HashSet<String> {
    map.iter()
        .filter(|(_, &v)| v >= thresh)
        .map(|(k, _)| k.clone())
        .collect()
}

/// Evaluate one sample. `present_thresh` is the minimum abundance to count as "detected".
pub fn evaluate(
    pred: &HashMap<String, f64>,
    truth: &HashMap<String, f64>,
    present_thresh: f64,
) -> Metrics {
    let p = presence(pred, present_thresh);
    let t = presence(truth, present_thresh);

    let tp = p.intersection(&t).count();
    let fp = p.difference(&t).count();
    let fn_ = t.difference(&p).count();

    let precision = if tp + fp == 0 {
        0.0
    } else {
        tp as f64 / (tp + fp) as f64
    };
    let recall = if tp + fn_ == 0 {
        0.0
    } else {
        tp as f64 / (tp + fn_) as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };

    // Abundance error over the union of labels.
    let mut labels: HashSet<&String> = HashSet::new();
    labels.extend(pred.keys());
    labels.extend(truth.keys());
    let mut l1 = 0.0;
    let mut num = 0.0; // sum |p-t|
    let mut den = 0.0; // sum (p+t)
    for l in labels {
        let pv = pred.get(l).copied().unwrap_or(0.0);
        let tv = truth.get(l).copied().unwrap_or(0.0);
        l1 += (pv - tv).abs();
        num += (pv - tv).abs();
        den += pv + tv;
    }
    let bray_curtis = if den == 0.0 { 0.0 } else { num / den };

    Metrics {
        tp,
        fp,
        fn_,
        precision,
        recall,
        f1,
        l1,
        bray_curtis,
    }
}

/// Parse a two-column `name<TAB>abundance` table (header lines starting with `#` ignored).
pub fn parse_abundance_tsv(text: &str) -> HashMap<String, f64> {
    let mut out = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(name), Some(val)) = (it.next(), it.next()) {
            if let Ok(v) = val.trim().parse::<f64>() {
                out.insert(name.trim().to_string(), v);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(pairs: &[(&str, f64)]) -> HashMap<String, f64> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    }

    #[test]
    fn perfect_prediction_scores_one() {
        let truth = m(&[("A", 0.7), ("B", 0.3)]);
        let pred = m(&[("A", 0.7), ("B", 0.3)]);
        let r = evaluate(&pred, &truth, 0.01);
        assert_eq!((r.precision, r.recall, r.f1), (1.0, 1.0, 1.0));
        assert!(r.l1 < 1e-9 && r.bray_curtis < 1e-9);
    }

    #[test]
    fn false_positive_and_abundance_error() {
        let truth = m(&[("A", 0.7), ("B", 0.3)]);
        let pred = m(&[("A", 0.6), ("B", 0.3), ("C", 0.1)]); // C is a false positive
        let r = evaluate(&pred, &truth, 0.05);
        assert_eq!((r.tp, r.fp, r.fn_), (2, 1, 0));
        assert!((r.recall - 1.0).abs() < 1e-9);
        assert!((r.precision - 2.0 / 3.0).abs() < 1e-9);
        // L1 = |0.6-0.7| + |0.3-0.3| + |0.1-0| = 0.2
        assert!((r.l1 - 0.2).abs() < 1e-9);
    }

    #[test]
    fn parse_tsv_skips_header() {
        let t = "#name\tabundance\nA\t0.7\nB\t0.3\n";
        let map = parse_abundance_tsv(t);
        assert_eq!(map.len(), 2);
        assert!((map["A"] - 0.7).abs() < 1e-9);
    }
}
