//! strain2bscan CLI.

mod cli;
mod commands;
mod panel;
mod report;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use cli::{Cli, Commands};
use strain2bscan::db::StrainDb;
use strain2bscan::enzymes::{parse_enzyme_set, Enzyme};
use strain2bscan::markers::{
    fastx_stem, genome_kmer_counts, genome_marker_counts_multi, is_fasta_path, kmer_db_token,
    parse_kmer_db_token, read_fastx, sample_kmer_counts_stream, sample_marker_counts_stream,
    single_copy_markers, sketch_threshold, Marker, MarkerCounts,
};
use strain2bscan::parallel::par_map;
use strain2bscan::quality::{self, GenomeRec, QualityFilter};

fn main() -> ExitCode {
    let cli = match Cli::try_parse() {
        Ok(c) => c,
        Err(e) => {
            let code = if e.use_stderr() { 2 } else { 0 };
            e.print().expect("failed to print clap error");
            return ExitCode::from(code);
        }
    };

    let result = match cli.command {
        Commands::Build(args) => commands::build::run(&args),
        Commands::Cluster(args) => commands::cluster::run(&args),
        Commands::DiagnoseTree(args) => commands::diagnose_tree::run(&args),
        Commands::Profile(args) => commands::profile::run(&args),
        Commands::MultiProfile(args) => commands::multi_profile::run(&args),
        Commands::Batch(args) => commands::batch::run(&args),
        Commands::Info(args) => commands::info::run(&args),
        Commands::Evaluate(args) => commands::evaluate::run(&args),
        Commands::Demo => commands::demo::run(),
        Commands::CstDemo => commands::cst_demo::run(),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Where markers come from: enzyme-digested 2bRAD tags (default), or sketched raw k-mers.
#[derive(Debug)]
pub(crate) enum MarkerSource {
    Enzyme(Vec<&'static Enzyme>),
    Kmer { k: usize, scale: u64 },
}

impl MarkerSource {
    /// Digest all contigs of one genome into marker copy numbers.
    fn genome_counts(&self, seqs: &[Vec<u8>]) -> MarkerCounts {
        match self {
            MarkerSource::Enzyme(set) => genome_marker_counts_multi(seqs, set),
            MarkerSource::Kmer { k, scale } => {
                genome_kmer_counts(seqs, *k, sketch_threshold(*scale))
            }
        }
    }

    /// Digest a sample's reads into marker counts.
    pub(crate) fn sample_counts(&self, reads: &Path) -> Result<MarkerCounts, String> {
        match self {
            MarkerSource::Enzyme(set) => sample_marker_counts_stream(reads, set),
            MarkerSource::Kmer { k, scale } => {
                sample_kmer_counts_stream(reads, *k, sketch_threshold(*scale))
            }
        }
        .map_err(|e| e.to_string())
    }

    /// Digest a possibly paired-end sample: R1 and R2 counted separately and summed.
    pub(crate) fn sample_counts_paired(
        &self,
        r1: &Path,
        r2: Option<&Path>,
    ) -> Result<MarkerCounts, String> {
        let mut counts = self.sample_counts(r1)?;
        if let Some(r2) = r2 {
            for (m, c) in self.sample_counts(r2)? {
                *counts.entry(m).or_insert(0) += c;
            }
        }
        Ok(counts)
    }

    /// What the DB header records in the `enzyme_csv` position.
    fn db_token(&self) -> Vec<String> {
        match self {
            MarkerSource::Enzyme(set) => enzyme_names(set),
            MarkerSource::Kmer { k, scale } => vec![kmer_db_token(*k, *scale)],
        }
    }

    /// Human-readable form for progress lines.
    pub(crate) fn describe(&self) -> String {
        match self {
            MarkerSource::Enzyme(set) => format!("enzymes: {}", enzyme_names(set).join("+")),
            MarkerSource::Kmer { k, scale } => format!("k-mer sketch: k={k}, scale={scale}"),
        }
    }

    /// Noun for the per-genome marker report ("tags" / "k-mers").
    fn marker_noun(&self) -> &'static str {
        match self {
            MarkerSource::Enzyme(_) => "tags",
            MarkerSource::Kmer { .. } => "k-mers",
        }
    }
}

fn enzyme_set(spec: &str) -> Result<Vec<&'static Enzyme>, String> {
    parse_enzyme_set(spec).ok_or_else(|| format!("unknown enzyme set: {spec}"))
}

fn enzyme_names(set: &[&Enzyme]) -> Vec<String> {
    set.iter().map(|e| e.name.to_string()).collect()
}

/// Parse `--kmer-size` (default 31) and `--sketch-scale` (default 100).
pub(crate) fn kmer_params(
    kmer_size: Option<usize>,
    sketch_scale: Option<u64>,
) -> Result<(usize, u64), String> {
    let k = kmer_size.unwrap_or(31);
    let scale = sketch_scale.unwrap_or(100);
    if k == 0 {
        return Err("bad --kmer-size (want integer >= 1)".into());
    }
    if scale == 0 {
        return Err("bad --sketch-scale (want integer >= 1)".into());
    }
    Ok((k, scale))
}

/// Marker source for the genome-side commands (`build`/`cluster`/`diagnose-tree`).
pub(crate) fn marker_source_for_build(
    marker_source: &str,
    enzyme: Option<&str>,
    kmer_size: usize,
    sketch_scale: u64,
) -> Result<MarkerSource, String> {
    match marker_source {
        "enzyme" => {
            let Some(spec) = enzyme else {
                return Err("missing --enzyme".into());
            };
            Ok(MarkerSource::Enzyme(enzyme_set(spec)?))
        }
        "kmer" => {
            if let Some(e) = enzyme {
                eprintln!("warning: --enzyme {e} is ignored with --marker-source kmer");
            }
            Ok(MarkerSource::Kmer {
                k: kmer_size,
                scale: sketch_scale,
            })
        }
        x => Err(format!("bad --marker-source {x} (want enzyme|kmer)")),
    }
}

/// Marker source for the sample side (`profile`): auto-detected from the DB header.
pub(crate) fn resolve_sample_source(
    db: &StrainDb,
    marker_source_arg: Option<&str>,
    enzyme_arg: Option<&str>,
    kmer_size_arg: Option<usize>,
    sketch_scale_arg: Option<u64>,
) -> Result<MarkerSource, String> {
    let db_kmer = if db.enzymes.len() == 1 {
        parse_kmer_db_token(&db.enzymes[0])
    } else {
        None
    };
    match (db_kmer, marker_source_arg) {
        (Some((k, scale)), None | Some("kmer")) => {
            if let Some(want_k) = kmer_size_arg {
                if want_k != k {
                    return Err(format!(
                        "--kmer-size {want_k} does not match the database (built with k={k})"
                    ));
                }
            }
            if let Some(want_scale) = sketch_scale_arg {
                if want_scale != scale {
                    return Err(format!(
                        "--sketch-scale {want_scale} does not match the database (built with scale={scale})"
                    ));
                }
            }
            if let Some(e) = enzyme_arg {
                eprintln!("warning: --enzyme {e} is ignored: the DB is a k-mer sketch (k={k}, scale={scale})");
            }
            Ok(MarkerSource::Kmer { k, scale })
        }
        (Some((k, scale)), Some("enzyme")) => Err(format!(
            "this database is a k-mer sketch (k={k}, scale={scale}); profiling it with \
             --marker-source enzyme would compare disjoint marker spaces"
        )),
        (Some(_), Some(x)) => Err(format!("bad --marker-source {x} (want enzyme|kmer)")),
        (None, Some("kmer")) => Err(format!(
            "this database was built from enzyme-digested tags ({}); --marker-source kmer \
             would compare disjoint marker spaces",
            db.enzymes.join("+")
        )),
        (None, None | Some("enzyme")) => {
            let set: Vec<&Enzyme> = if !db.enzymes.is_empty() {
                parse_enzyme_set(&db.enzymes.join(",")).ok_or("DB records an unknown enzyme")?
            } else {
                let Some(spec) = enzyme_arg else {
                    return Err("missing --enzyme".into());
                };
                enzyme_set(spec)?
            };
            Ok(MarkerSource::Enzyme(set))
        }
        (None, Some(x)) => Err(format!("bad --marker-source {x} (want enzyme|kmer)")),
    }
}

/// Digest every FASTA genome in `dir` → `GenomeRec`, in parallel across genomes.
pub(crate) fn digest_genome_dir(
    dir: &Path,
    source: &MarkerSource,
) -> Result<Vec<GenomeRec>, String> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|e| e.to_string())? {
        let path = entry.map_err(|e| e.to_string())?.path();
        if path.is_file() && is_fasta_path(&path) {
            paths.push(path);
        }
    }
    if paths.is_empty() {
        return Err("no FASTA genomes (.fa/.fasta/.fna, optionally .gz) found".into());
    }
    paths.sort();
    let results: Vec<Result<GenomeRec, String>> = par_map(&paths, |path| {
        let name = fastx_stem(path);
        let seqs = read_fastx(path).map_err(|e| e.to_string())?;
        let n_contigs = seqs.len();
        let counts = source.genome_counts(&seqs);
        let full_markers: Vec<Marker> = counts.keys().copied().collect();
        Ok(GenomeRec {
            name,
            n_contigs,
            markers: single_copy_markers(&counts),
            full_markers,
        })
    });
    results.into_iter().collect()
}

/// Parse `--max-contigs` / `--min-tag-fraction` into a `QualityFilter`.
pub(crate) fn parse_quality_filter(
    max_contigs: Option<usize>,
    min_tag_fraction: Option<f64>,
) -> Result<QualityFilter, String> {
    Ok(QualityFilter {
        max_contigs,
        min_tag_fraction,
        ..QualityFilter::default()
    })
}

/// Digest a genome dir, apply the assembly-quality filter, print the report, and return the kept recs.
pub(crate) fn digest_and_filter(
    dir: &Path,
    source: &MarkerSource,
    max_contigs: Option<usize>,
    min_tag_fraction: Option<f64>,
) -> Result<Vec<GenomeRec>, String> {
    let genomes = digest_genome_dir(dir, source)?;
    let filt = parse_quality_filter(max_contigs, min_tag_fraction)?;
    let rep = quality::apply(genomes, &filt);
    println!(
        "quality: {} genomes, median single-copy {} = {}",
        rep.n_input,
        source.marker_noun(),
        rep.median_tags
    );
    for (name, nt) in &rep.flagged {
        println!(
            "  ⚠ likely incomplete: {name} has {nt} tags (< {:.0}% of median {}) — kept; pass --min-tag-fraction to drop",
            filt.warn_fraction * 100.0,
            rep.median_tags
        );
    }
    for (name, reason) in &rep.dropped {
        println!("  ✗ dropped {name}: {reason}");
    }
    if rep.kept.is_empty() {
        return Err("all genomes removed by the quality filter".into());
    }
    Ok(rep.kept)
}

pub(crate) fn print_stats(db: &StrainDb) {
    let s = db.stats();
    println!(
        "  units={}  markers={}  unique={} ({:.1}%)  avg_markers/unit={:.0}",
        s.n_strains,
        s.n_markers,
        s.unique_markers,
        s.unique_fraction * 100.0,
        s.avg_markers_per_strain
    );
}

#[cfg(test)]
mod tests {
    use super::{resolve_sample_source, MarkerSource};
    use strain2bscan::db::StrainDb;

    fn kmer_db() -> StrainDb {
        let mut db = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
        ]);
        db.enzymes = vec!["kmer15s1".to_string()];
        db
    }

    fn enzyme_db() -> StrainDb {
        let mut db = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
        ]);
        db.enzymes = vec!["BcgI".to_string()];
        db
    }

    #[test]
    fn profile_auto_detects_kmer_db() {
        let src = resolve_sample_source(&kmer_db(), None, None, None, None).unwrap();
        match src {
            MarkerSource::Kmer { k, scale } => assert_eq!((k, scale), (15, 1)),
            _ => panic!("k-mer DB must select the k-mer source"),
        }
        let src = resolve_sample_source(&kmer_db(), Some("kmer"), None, Some(15), Some(1)).unwrap();
        assert!(matches!(src, MarkerSource::Kmer { k: 15, scale: 1 }));
    }

    #[test]
    fn profile_rejects_kmer_param_mismatch() {
        let err = resolve_sample_source(&kmer_db(), None, None, Some(31), None).unwrap_err();
        assert!(err.contains("does not match the database"), "{err}");
        let err = resolve_sample_source(&kmer_db(), None, None, None, Some(100)).unwrap_err();
        assert!(err.contains("does not match the database"), "{err}");
    }

    #[test]
    fn profile_rejects_crossing_marker_spaces() {
        let err = resolve_sample_source(&kmer_db(), Some("enzyme"), None, None, None).unwrap_err();
        assert!(err.contains("k-mer sketch"), "{err}");
        let err = resolve_sample_source(&enzyme_db(), Some("kmer"), None, None, None).unwrap_err();
        assert!(err.contains("disjoint marker spaces"), "{err}");
        let err = resolve_sample_source(&kmer_db(), Some("wat"), None, None, None).unwrap_err();
        assert!(err.contains("bad --marker-source"), "{err}");
    }

    #[test]
    fn profile_enzyme_db_stays_on_enzyme_path() {
        let src = resolve_sample_source(&enzyme_db(), None, None, None, None).unwrap();
        match src {
            MarkerSource::Enzyme(set) => {
                assert_eq!(set.len(), 1);
                assert_eq!(set[0].name, "BcgI");
            }
            _ => panic!("enzyme DB must select the enzyme source"),
        }
    }
}
