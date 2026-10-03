//! Derive-based clap CLI for strain2bscan.

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "strain2bscan",
    version,
    about = "Fast strain-level metagenomic profiling on 2bRAD-reduced k-mer markers",
    long_about = "Strain2bScan: Fast strain-level metagenomic profiling on 2bRAD-reduced k-mer markers (a 2bRAD extension of StrainScan).\n\
        \n\
        `<set>` is `all` (all 16 type-IIB enzymes), a single enzyme (`BcgI`), or a comma list \
        (`BcgI,CspCI`). Use `BcgI` for BcgI 2bRAD data; use `all` to digitally digest a \
        conventional metagenome and enrich strain-specific markers. The genome DB and the sample \
        must use the same enzyme set — `profile` reads the set from the DB header automatically.\n\
        \n\
        `--marker-source kmer` replaces enzyme digestion with sketched canonical k-mers \
        (k = `--kmer-size`, default 31; keep ~1/`--sketch-scale`, default 100, of them by hash). \
        A k-mer DB records `kmer<K>s<S>` in its header; `profile`/`multi-profile` auto-detect \
        this and refuse to cross the two marker spaces."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Build a strain DB from genomes in a directory.
    Build(BuildArgs),
    /// Build a within-species Cluster Search Tree DB from genomes.
    Cluster(ClusterArgs),
    /// Measure whether a hierarchical Cluster Search Tree could work on this panel.
    DiagnoseTree(DiagnoseTreeArgs),
    /// Profile a single species sample against a strain DB.
    Profile(ProfileArgs),
    /// Multi-species strain profiling against per-species DBs.
    #[command(name = "multi-profile")]
    MultiProfile(MultiProfileArgs),
    /// Profile many samples against a multi-species panel from a manifest.
    Batch(BatchArgs),
    /// Print information about a strain DB.
    Info(InfoArgs),
    /// Evaluate predictions against a truth table.
    Evaluate(EvaluateArgs),
    /// In-memory conspecific demo.
    Demo,
    /// Cluster Search Tree demo.
    #[command(name = "cst-demo")]
    CstDemo,
}

#[derive(Parser)]
pub struct BuildArgs {
    /// Directory containing FASTA genomes (.fa/.fasta/.fna, optionally .gz).
    #[arg(long)]
    pub genomes: PathBuf,
    /// Enzyme set: `all`, `recommended`, a single enzyme (`BcgI`), or a comma list.
    #[arg(long)]
    pub enzyme: Option<String>,
    /// Output DB path.
    #[arg(long)]
    pub out: PathBuf,
    /// Max contigs for assembly quality filter.
    #[arg(long)]
    pub max_contigs: Option<usize>,
    /// Min fraction of median tags for assembly quality filter.
    #[arg(long)]
    pub min_tag_fraction: Option<f64>,
    /// Marker source: `enzyme` (default) or `kmer`.
    #[arg(long, default_value = "enzyme")]
    pub marker_source: String,
    /// K-mer size (only used with `--marker-source kmer`).
    #[arg(long, default_value_t = 31)]
    pub kmer_size: usize,
    /// Sketch scale denominator (only used with `--marker-source kmer`).
    #[arg(long, default_value_t = 100)]
    pub sketch_scale: u64,
}

#[derive(Parser)]
pub struct ClusterArgs {
    #[arg(long)]
    pub genomes: PathBuf,
    #[arg(long)]
    pub enzyme: Option<String>,
    #[arg(long)]
    pub out: PathBuf,
    /// Clustering similarity threshold.
    #[arg(long, default_value_t = 0.95)]
    pub similarity: f64,
    /// Use max-containment instead of Jaccard for uneven-completeness panels.
    #[arg(long)]
    pub containment: bool,
    #[arg(long)]
    pub max_contigs: Option<usize>,
    #[arg(long)]
    pub min_tag_fraction: Option<f64>,
    #[arg(long, default_value = "enzyme")]
    pub marker_source: String,
    #[arg(long, default_value_t = 31)]
    pub kmer_size: usize,
    #[arg(long, default_value_t = 100)]
    pub sketch_scale: u64,
}

#[derive(Parser)]
pub struct DiagnoseTreeArgs {
    #[arg(long)]
    pub genomes: PathBuf,
    #[arg(long)]
    pub enzyme: Option<String>,
    #[arg(long, default_value_t = 0.95)]
    pub similarity: f64,
    #[arg(long)]
    pub containment: bool,
    #[arg(long)]
    pub max_contigs: Option<usize>,
    #[arg(long)]
    pub min_tag_fraction: Option<f64>,
    #[arg(long, default_value = "enzyme")]
    pub marker_source: String,
    #[arg(long, default_value_t = 31)]
    pub kmer_size: usize,
    #[arg(long, default_value_t = 100)]
    pub sketch_scale: u64,
}

#[derive(Parser)]
pub struct ProfileArgs {
    /// Strain DB path.
    #[arg(long)]
    pub db: PathBuf,
    /// Sample reads (.fa/.fasta/.fna/.fq/.fastq, optionally .gz).
    #[arg(long)]
    pub reads: PathBuf,
    /// Enzyme set (auto-detected from DB if omitted; required for enzyme DB if not recorded).
    #[arg(long)]
    pub enzyme: Option<String>,
    /// Output predictions TSV.
    #[arg(long)]
    pub out: Option<PathBuf>,
    /// Min supporting unique markers (default 8).
    #[arg(long)]
    pub min_support: Option<usize>,
    /// Min fraction of unique markers detected (default 0.1).
    #[arg(long)]
    pub min_coverage: Option<f64>,
    /// Min within-species relative abundance (default 0.0).
    #[arg(long)]
    pub min_abundance: Option<f64>,
    /// Min abundance ratio between consecutive calls to cut trace tail (default 0.0 = off).
    #[arg(long)]
    pub trace_gap: Option<f64>,
    /// Absolute abundance floor applied after trace-gap cut (default 0.0).
    #[arg(long)]
    pub trace_floor: Option<f64>,
    /// Min coverage / (1 - e^(-depth)) consistency (default 0.5).
    #[arg(long)]
    pub min_consistency: Option<f64>,
    /// ElasticNet alpha (default 0.0).
    #[arg(long)]
    pub enet_alpha: Option<f64>,
    /// Layer-1 algorithm: `auto`, `unique`, or `cst` (default auto).
    #[arg(long)]
    pub layer1: Option<String>,
    /// Layer-2 algorithm: `depth` or `enet` (default depth).
    #[arg(long)]
    pub layer2: Option<String>,
    /// Restore pre-0.2 fixed singleton filter and unscaled species floor.
    #[arg(long)]
    pub fixed_gate: bool,
    /// Marker source override: `enzyme` or `kmer`.
    #[arg(long)]
    pub marker_source: Option<String>,
    #[arg(long)]
    pub kmer_size: Option<usize>,
    #[arg(long)]
    pub sketch_scale: Option<u64>,
}

#[derive(Parser)]
pub struct MultiProfileArgs {
    /// Directory containing per-species DBs (*.tsv).
    #[arg(long)]
    pub dbs: PathBuf,
    /// Sample reads.
    #[arg(long)]
    pub reads: PathBuf,
    /// Enzyme set (required for enzyme panels; ignored for k-mer panels).
    #[arg(long)]
    pub enzyme: Option<String>,
    #[arg(long)]
    pub out: Option<PathBuf>,
    /// Min species-specific markers to attempt strain resolution (default 200).
    #[arg(long)]
    pub min_species_markers: Option<usize>,
    /// Min fraction of species-specific panel required (default 0.0).
    #[arg(long)]
    pub min_species_marker_frac: Option<f64>,
    /// Min detected species-specific markers to call presence (default 10).
    #[arg(long)]
    pub min_species_detect: Option<usize>,
    /// Min within-species relative abundance.
    #[arg(long)]
    pub min_abundance: Option<f64>,
    /// Min cross-species global abundance to keep a call (default 0.0).
    #[arg(long)]
    pub min_global_abundance: Option<f64>,
    /// Min supporting unique markers (default 8).
    #[arg(long)]
    pub min_support: Option<usize>,
    /// Min fraction of unique markers detected (default 0.1).
    #[arg(long)]
    pub min_coverage: Option<f64>,
    #[arg(long)]
    pub trace_gap: Option<f64>,
    #[arg(long)]
    pub trace_floor: Option<f64>,
    #[arg(long)]
    pub min_consistency: Option<f64>,
    /// ElasticNet alpha (default 0.0).
    #[arg(long)]
    pub enet_alpha: Option<f64>,
    /// Layer-1 algorithm: `auto`, `unique`, or `cst` (default auto).
    #[arg(long)]
    pub layer1: Option<String>,
    /// Layer-2 algorithm: `depth` or `enet` (default depth).
    #[arg(long)]
    pub layer2: Option<String>,
    #[arg(long)]
    pub fixed_gate: bool,
    /// Disable adaptive singleton admission.
    #[arg(long)]
    pub no_adaptive_singleton: bool,
    /// Disable adaptive floor relaxation.
    #[arg(long)]
    pub no_adaptive_floor: bool,
    /// Disable cross-species marker filtering.
    #[arg(long)]
    pub no_cross_species_filter: bool,
    #[arg(long)]
    pub marker_source: Option<String>,
    #[arg(long)]
    pub kmer_size: Option<usize>,
    #[arg(long)]
    pub sketch_scale: Option<u64>,
}

#[derive(Parser)]
pub struct BatchArgs {
    #[arg(long)]
    pub dbs: PathBuf,
    /// Manifest CSV with header `sample,reads1[,reads2]`.
    #[arg(long)]
    pub manifest: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long)]
    pub enzyme: Option<String>,
    #[arg(long)]
    pub min_species_markers: Option<usize>,
    #[arg(long)]
    pub min_species_marker_frac: Option<f64>,
    #[arg(long)]
    pub min_species_detect: Option<usize>,
    #[arg(long)]
    pub min_abundance: Option<f64>,
    #[arg(long)]
    pub min_global_abundance: Option<f64>,
    /// Min supporting unique markers (default 8).
    #[arg(long)]
    pub min_support: Option<usize>,
    /// Min fraction of unique markers detected (default 0.1).
    #[arg(long)]
    pub min_coverage: Option<f64>,
    #[arg(long)]
    pub trace_gap: Option<f64>,
    #[arg(long)]
    pub trace_floor: Option<f64>,
    #[arg(long)]
    pub min_consistency: Option<f64>,
    /// ElasticNet alpha (default 0.0).
    #[arg(long)]
    pub enet_alpha: Option<f64>,
    /// Layer-1 algorithm: `auto`, `unique`, or `cst` (default auto).
    #[arg(long)]
    pub layer1: Option<String>,
    /// Layer-2 algorithm: `depth` or `enet` (default depth).
    #[arg(long)]
    pub layer2: Option<String>,
    #[arg(long)]
    pub fixed_gate: bool,
    #[arg(long)]
    pub no_adaptive_singleton: bool,
    #[arg(long)]
    pub no_adaptive_floor: bool,
    #[arg(long)]
    pub no_cross_species_filter: bool,
    #[arg(long)]
    pub marker_source: Option<String>,
    #[arg(long)]
    pub kmer_size: Option<usize>,
    #[arg(long)]
    pub sketch_scale: Option<u64>,
}

#[derive(Parser)]
pub struct InfoArgs {
    #[arg(long)]
    pub db: PathBuf,
}

#[derive(Parser)]
pub struct EvaluateArgs {
    #[arg(long)]
    pub pred: PathBuf,
    #[arg(long)]
    pub truth: PathBuf,
    /// Abundance threshold to count a strain as present.
    #[arg(long, default_value_t = 0.01)]
    pub present: f64,
}
