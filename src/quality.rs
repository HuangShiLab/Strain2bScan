//! Assembly-quality filtering for reference genomes.
//!
//! Variable assembly completeness biases Jaccard clustering toward **spurious splits**: an
//! incomplete genome's marker set is ~a subset of its complete twin's, so their Jaccard
//! (|A∩B|/|A∪B|) drops well below 1 and they fail to cluster even though they are the same
//! strain. We don't run CheckM inline (heavy external dependency); instead we use two
//! dependency-free proxies computed from data we already have:
//!   * **contig count** — fragmentation / possible contamination (`--max-contigs`)
//!   * **single-copy tag count vs the conspecific median** — completeness proxy
//!     (`--min-tag-fraction`): an incomplete genome yields proportionally fewer tags.
//!
//! Flagging (warn) is always on; dropping happens only when a threshold is set. This assumes
//! most genomes in the input are reasonably complete (the median is the reference).

use crate::markers::Marker;

#[derive(Debug, Clone)]
pub struct GenomeRec {
    pub name: String,
    pub n_contigs: usize,
    pub markers: Vec<Marker>,
}

#[derive(Debug, Clone)]
pub struct QualityFilter {
    /// Drop genomes with more than this many contigs (None = no contig filter).
    pub max_contigs: Option<usize>,
    /// Drop genomes whose tag count < fraction × median tag count (None = no completeness drop).
    pub min_tag_fraction: Option<f64>,
    /// Flag (warn, but keep) genomes whose tag count < this × median. Default 0.5.
    pub warn_fraction: f64,
}

impl Default for QualityFilter {
    fn default() -> Self {
        QualityFilter { max_contigs: None, min_tag_fraction: None, warn_fraction: 0.5 }
    }
}

#[derive(Debug)]
pub struct QualityReport {
    pub n_input: usize,
    pub median_tags: usize,
    /// kept-but-suspect (name, n_tags) — likely incomplete.
    pub flagged: Vec<(String, usize)>,
    /// removed (name, reason).
    pub dropped: Vec<(String, String)>,
    pub kept: Vec<GenomeRec>,
}

/// Median single-copy tag count across genomes (the completeness reference).
pub fn median_tags(genomes: &[GenomeRec]) -> usize {
    if genomes.is_empty() {
        return 0;
    }
    let mut v: Vec<usize> = genomes.iter().map(|g| g.markers.len()).collect();
    v.sort_unstable();
    v[v.len() / 2]
}

/// Apply the quality filter, returning kept genomes plus a flag/drop report.
pub fn apply(genomes: Vec<GenomeRec>, f: &QualityFilter) -> QualityReport {
    let n_input = genomes.len();
    let median = median_tags(&genomes);
    let drop_below = f.min_tag_fraction.map(|fr| (fr * median as f64) as usize);
    let warn_below = (f.warn_fraction * median as f64) as usize;

    let mut kept = Vec::new();
    let mut flagged = Vec::new();
    let mut dropped = Vec::new();
    for g in genomes {
        let nt = g.markers.len();
        let mut reason: Option<String> = None;
        if let Some(mc) = f.max_contigs {
            if g.n_contigs > mc {
                reason = Some(format!("{} contigs > --max-contigs {mc}", g.n_contigs));
            }
        }
        if reason.is_none() {
            if let (Some(db), Some(fr)) = (drop_below, f.min_tag_fraction) {
                if nt < db {
                    reason = Some(format!("{nt} tags < {:.0}% of median {median}", fr * 100.0));
                }
            }
        }
        match reason {
            Some(r) => dropped.push((g.name.clone(), r)),
            None => {
                if median > 0 && nt < warn_below {
                    flagged.push((g.name.clone(), nt));
                }
                kept.push(g);
            }
        }
    }
    QualityReport { n_input, median_tags: median, flagged, dropped, kept }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, n_contigs: usize, n_tags: usize) -> GenomeRec {
        GenomeRec { name: name.into(), n_contigs, markers: (0..n_tags as u64).collect() }
    }

    #[test]
    fn median_is_robust_to_a_few_incompletes() {
        let g = vec![rec("a", 1, 1000), rec("b", 1, 1010), rec("c", 1, 990), rec("d", 50, 300)];
        assert_eq!(median_tags(&g), 1000);
    }

    #[test]
    fn flags_low_completeness_without_dropping_by_default() {
        let g = vec![rec("complete1", 1, 1000), rec("complete2", 1, 1000), rec("partial", 1, 300)];
        let r = apply(g, &QualityFilter::default());
        assert_eq!(r.kept.len(), 3); // nothing dropped by default
        assert_eq!(r.dropped.len(), 0);
        assert_eq!(r.flagged.len(), 1); // partial (300 < 0.5*1000) flagged
        assert_eq!(r.flagged[0].0, "partial");
    }

    #[test]
    fn drops_on_min_tag_fraction_and_max_contigs() {
        let g = vec![
            rec("good", 2, 1000),
            rec("good2", 3, 1000),
            rec("incomplete", 4, 200),  // 200 < 0.5*1000 -> dropped
            rec("fragmented", 800, 1000), // contigs > 500 -> dropped
        ];
        let f = QualityFilter { max_contigs: Some(500), min_tag_fraction: Some(0.5), warn_fraction: 0.5 };
        let r = apply(g, &f);
        let kept: Vec<&str> = r.kept.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(kept, vec!["good", "good2"]);
        assert_eq!(r.dropped.len(), 2);
    }
}
