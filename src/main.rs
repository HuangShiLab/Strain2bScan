//! strain2bscan CLI (prototype).
//!
//!   strain2bscan build    --genomes <dir> --enzyme <set> --out <db.tsv> [--max-contigs N] [--min-tag-fraction F]
//!   strain2bscan cluster  --genomes <dir> --enzyme <set> --out <clusterdb.tsv> [--similarity 0.95] [--containment (uneven-completeness panels)] [--max-contigs N] [--min-tag-fraction F]
//!   strain2bscan profile  --db <db.tsv> --reads <fastx> [--enzyme <set>] [--out pred.tsv] [--min-support N] [--min-coverage F] [--min-abundance F] [--fixed-gate]
//!   strain2bscan multi-profile --dbs <dir> --reads <fastx> --enzyme <set> [--out pred.tsv] [--fixed-gate]
//!   strain2bscan info     --db <db.tsv>
//!   strain2bscan evaluate --pred <pred.tsv> --truth <truth.tsv> [--present 0.01]
//!   strain2bscan demo | cst-demo
//!
//! `<set>` is `all` (all 16 type-IIB enzymes), a single enzyme (`BcgI`), or a comma list
//! (`BcgI,CspCI`). Use `BcgI` for BcgI 2bRAD data; use `all` to digitally digest a
//! conventional metagenome and enrich strain-specific markers. The genome DB and the sample
//! must use the same enzyme set — `profile` reads the set from the DB header automatically.
//!
//! Arg parsing is hand-rolled to keep the prototype dependency-free; production uses clap.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use strain2bscan::bench::{evaluate, parse_abundance_tsv};
use strain2bscan::cst::{SpeciesCst, DEFAULT_SIMILARITY, MINHASH_ABOVE};
use strain2bscan::db::StrainDb;
use strain2bscan::enzymes::{parse_enzyme_set, Enzyme};
use strain2bscan::fxhash::{FxHashMap, FxHashSet};
use strain2bscan::identify::{detectable_fraction, naive_profile, profile, Params, StrainCall};
use strain2bscan::markers::{
    fastx_stem, genome_marker_counts_multi, is_fasta_path, read_fastx, sample_marker_counts_stream,
    single_copy_markers, Marker, MarkerCounts,
};
use strain2bscan::parallel::{num_threads, par_map};
use strain2bscan::quality::{self, GenomeRec, QualityFilter};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("");
    let opts = parse_opts(&args);

    let result = match cmd {
        "build" => cmd_build(&opts),
        "cluster" => cmd_cluster(&opts),
        "profile" => cmd_profile(&opts),
        "multi-profile" => cmd_multi_profile(&opts),
        "info" => cmd_info(&opts),
        "evaluate" => cmd_evaluate(&opts),
        "demo" => cmd_demo(),
        "cst-demo" => cmd_cst_demo(),
        _ => {
            eprintln!(
                "usage:\n  \
                 strain2bscan build    --genomes <dir> --enzyme <set> --out <db.tsv> [--max-contigs N] [--min-tag-fraction F]\n  \
                 strain2bscan cluster  --genomes <dir> --enzyme <set> --out <clusterdb.tsv> [--similarity 0.95] [--containment (uneven-completeness panels)] [--max-contigs N] [--min-tag-fraction F]\n  \
                 strain2bscan profile  --db <db.tsv> --reads <fastx> [--enzyme <set>] [--out pred.tsv] [--min-support N] [--min-coverage F] [--min-abundance F] [--fixed-gate]\n  \
                 strain2bscan multi-profile --dbs <dir> --reads <fastx> --enzyme <set> [--out pred.tsv] [--min-species-markers N] [--min-species-marker-frac F] [--min-species-detect N] [--min-abundance F] [--fixed-gate|--no-adaptive-singleton|--no-adaptive-floor] [--no-cross-species-filter]   (many species, sample digested once)\n  \
                 strain2bscan info     --db <db.tsv>\n  \
                 strain2bscan evaluate --pred <pred.tsv> --truth <truth.tsv> [--present 0.01]\n  \
                 strain2bscan demo | cst-demo\n\n\
                 <set> = all | BcgI | BcgI,CspCI  (use BcgI for 2bRAD data; all for conventional metagenomes)\n\
                 reads/genomes may be gzipped (.fq.gz, .fna.gz).\n\
                 multi-profile reports two scopes: `abundance` sums to 1.0 WITHIN each species\n\
                 (the primary number), `global_abundance` sums to 1.0 over the strain-resolved\n\
                 part of the sample. --min-abundance applies within-species. Detection and depth\n\
                 use only markers specific to their species across the whole panel, so a\n\
                 co-present congener cannot inflate a cluster's depth; disable with\n\
                 --no-cross-species-filter (single-species `profile` cannot do this).\n\
                 --fixed-gate restores the pre-0.2 fixed singleton filter (count>=2) AND the\n\
                 unscaled species floor. The two can be reverted independently with\n\
                 --no-adaptive-singleton / --no-adaptive-floor: on a large panel the floor\n\
                 relaxation costs specificity (an ABSENT species has the lowest estimated depth,\n\
                 so it receives the largest relaxation), while the singleton admission is what\n\
                 buys low-depth sensitivity. Calibrate the pair on your own panel."
            );
            return ExitCode::from(2);
        }
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Collect `--key value` pairs (and bare `--flag`).
fn parse_opts(args: &[String]) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let mut i = 1;
    while i < args.len() {
        if let Some(key) = args[i].strip_prefix("--") {
            let val = args.get(i + 1).filter(|v| !v.starts_with("--"));
            match val {
                Some(v) => {
                    map.insert(key.to_string(), v.clone());
                    i += 2;
                }
                None => {
                    map.insert(key.to_string(), "true".into());
                    i += 1;
                }
            }
        } else {
            i += 1;
        }
    }
    map
}

fn req<'a>(opts: &'a HashMap<String, String>, key: &str) -> Result<&'a String, String> {
    opts.get(key).ok_or_else(|| format!("missing --{key}"))
}

fn enzyme_set(opts: &HashMap<String, String>) -> Result<Vec<&'static Enzyme>, String> {
    let spec = req(opts, "enzyme")?;
    parse_enzyme_set(spec).ok_or_else(|| format!("unknown enzyme set: {spec}"))
}

fn enzyme_names(set: &[&Enzyme]) -> Vec<String> {
    set.iter().map(|e| e.name.to_string()).collect()
}

/// Digest every FASTA genome in `dir` → `GenomeRec` (name, contig count, single-copy tag
/// markers), in parallel across genomes (the dominant build cost).
fn digest_genome_dir(dir: &Path, enzymes: &[&Enzyme]) -> Result<Vec<GenomeRec>, String> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|e| e.to_string())? {
        let path = entry.map_err(|e| e.to_string())?.path();
        // `.fna.gz` included: RefSeq/ENA assemblies ship gzipped, and decompressing a whole
        // reference panel just to build a DB is pure friction. FASTA only — a stray read file
        // in this directory must not be digested as though it were a reference genome.
        if path.is_file() && is_fasta_path(&path) {
            paths.push(path);
        }
    }
    if paths.is_empty() {
        return Err("no FASTA genomes (.fa/.fasta/.fna, optionally .gz) found".into());
    }
    paths.sort(); // deterministic genome order regardless of threading
    let results: Vec<Result<GenomeRec, String>> = par_map(&paths, |path| {
        // `fastx_stem`, not `file_stem`: the latter leaves `X.fna.gz` named `X.fna`, so the same
        // genome would get a different identifier depending only on whether it was compressed.
        let name = fastx_stem(path);
        let seqs = read_fastx(path).map_err(|e| e.to_string())?;
        let n_contigs = seqs.len();
        let counts = genome_marker_counts_multi(&seqs, enzymes);
        let full_markers: Vec<Marker> = counts.keys().copied().collect();
        Ok(GenomeRec { name, n_contigs, markers: single_copy_markers(&counts), full_markers })
    });
    results.into_iter().collect()
}

/// Parse `--max-contigs` / `--min-tag-fraction` into a `QualityFilter`.
fn parse_quality_filter(opts: &HashMap<String, String>) -> Result<QualityFilter, String> {
    let max_contigs = match opts.get("max-contigs") {
        Some(s) => Some(s.parse().map_err(|_| "bad --max-contigs (want integer)")?),
        None => None,
    };
    let min_tag_fraction = match opts.get("min-tag-fraction") {
        Some(s) => Some(s.parse().map_err(|_| "bad --min-tag-fraction (want 0..1)")?),
        None => None,
    };
    Ok(QualityFilter { max_contigs, min_tag_fraction, ..QualityFilter::default() })
}

/// Digest a genome dir, apply the assembly-quality filter, print the report, and return the
/// kept `(name, markers)`. Variable assembly completeness biases Jaccard clustering toward
/// spurious splits; flagging is always on, dropping happens only when a threshold is set.
fn digest_and_filter(
    dir: &Path,
    enzymes: &[&Enzyme],
    opts: &HashMap<String, String>,
) -> Result<Vec<GenomeRec>, String> {
    let genomes = digest_genome_dir(dir, enzymes)?;
    let filt = parse_quality_filter(opts)?;
    let rep = quality::apply(genomes, &filt);
    println!(
        "quality: {} genomes, median single-copy tags = {}",
        rep.n_input, rep.median_tags
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

fn cmd_build(opts: &HashMap<String, String>) -> Result<(), String> {
    let set = enzyme_set(opts)?;
    let genomes = PathBuf::from(req(opts, "genomes")?);
    let out = PathBuf::from(req(opts, "out")?);

    let recs = digest_and_filter(&genomes, &set, opts)?;
    for r in &recs {
        println!("  {}: {} single-copy tag markers", r.name, r.markers.len());
    }
    let mut db = StrainDb::build(recs.into_iter().map(|r| (r.name, r.markers)).collect());
    db.enzymes = enzyme_names(&set);
    db.save(&out).map_err(|e| e.to_string())?;
    print_stats(&db);
    println!("saved DB ({}) -> {}", db.enzymes.join("+"), out.display());
    Ok(())
}

/// Build a within-species Cluster Search Tree DB from genomes (StrainScan Layer-1/2 step).
fn cmd_cluster(opts: &HashMap<String, String>) -> Result<(), String> {
    let set = enzyme_set(opts)?;
    let genomes = PathBuf::from(req(opts, "genomes")?);
    let out = PathBuf::from(req(opts, "out")?);
    let similarity = opts
        .get("similarity")
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_SIMILARITY);

    let containment = opts.contains_key("containment");
    let recs = digest_and_filter(&genomes, &set, opts)?;
    let n_genomes = recs.len();
    let cst = SpeciesCst::build(
        recs.into_iter().map(|r| (r.name, r.markers, r.full_markers)).collect(),
        similarity,
        containment,
    );
    let dist = if containment {
        "max-containment"
    } else {
        "Jaccard"
    };
    let method = if n_genomes > MINHASH_ABOVE {
        "MinHash"
    } else {
        "exact"
    };
    println!(
        "clustered {} genomes into {} cluster(s) @ similarity {similarity} (enzymes: {}, threads: {}, clustering: {}-{})",
        cst.genome_names.len(),
        cst.n_clusters(),
        enzyme_names(&set).join("+"),
        num_threads(),
        method,
        dist
    );
    for (cid, members) in cst.clusters.iter().enumerate() {
        let names: Vec<&str> = members
            .iter()
            .map(|&g| cst.genome_names[g].as_str())
            .collect();
        println!("  C{cid}: {}", names.join(", "));
    }
    let s = cst.marker_class_summary();
    println!(
        "  marker classes: species_core={}  shared_partial={}  cluster_specific={}  strain_specific={}",
        s.get("species_core").unwrap_or(&0),
        s.get("shared_partial").unwrap_or(&0),
        s.get("cluster_specific").unwrap_or(&0),
        s.get("strain_specific").unwrap_or(&0),
    );

    // Resolvability check: a cluster needs enough cluster-specific markers to be detectable.
    let mut db = cst.cluster_db();
    db.enzymes = enzyme_names(&set);
    let min_markers = Params::default().min_support_markers;
    let mut resolvable = 0usize;
    for cid in 0..db.n_strains() {
        let n_spec = db.unique_marker_count(cid);
        if n_spec >= min_markers {
            resolvable += 1;
        } else {
            println!(
                "  ⚠ C{cid} has only {n_spec} cluster-specific markers (< {min_markers}); \
                 not reliably resolvable with this enzyme set."
            );
        }
    }
    if resolvable == 0 {
        println!(
            "  ✗ NOT DOABLE at strain/cluster level for this species with enzyme(s) {}. \
             The species can still be detected (Layer-1); for finer resolution use more \
             enzymes (--enzyme all) on a conventional metagenome.",
            enzyme_names(&set).join("+")
        );
    }

    // Write membership sidecar (genome -> cluster) for benchmark truth remapping.
    let members_path = out.with_extension("members.tsv");
    {
        use std::io::Write;
        let mut w = std::fs::File::create(&members_path).map_err(|e| e.to_string())?;
        writeln!(w, "#genome\tcluster").map_err(|e| e.to_string())?;
        for (cid, members) in cst.clusters.iter().enumerate() {
            for &g in members {
                writeln!(w, "{}\tC{cid}", cst.genome_names[g]).map_err(|e| e.to_string())?;
            }
        }
    }

    db.save(&out).map_err(|e| e.to_string())?;
    println!(
        "saved cluster DB -> {} ({} clusters, {} resolvable); membership -> {}",
        out.display(),
        db.n_strains(),
        resolvable,
        members_path.display()
    );
    Ok(())
}

fn cmd_evaluate(opts: &HashMap<String, String>) -> Result<(), String> {
    let pred_text = std::fs::read_to_string(req(opts, "pred")?).map_err(|e| e.to_string())?;
    let truth_text = std::fs::read_to_string(req(opts, "truth")?).map_err(|e| e.to_string())?;
    let present = opts
        .get("present")
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.01);
    let pred = parse_abundance_tsv(&pred_text);
    let truth = parse_abundance_tsv(&truth_text);
    let m = evaluate(&pred, &truth, present);
    println!(
        "TP={} FP={} FN={}  precision={:.3} recall={:.3} F1={:.3}  L1={:.3} Bray-Curtis={:.3}",
        m.tp, m.fp, m.fn_, m.precision, m.recall, m.f1, m.l1, m.bray_curtis
    );
    Ok(())
}

/// Parse the Layer-2 tuning flags shared by `profile` and `multi-profile`.
fn parse_params(opts: &HashMap<String, String>) -> Result<Params, String> {
    let mut p = Params::default();
    if let Some(v) = opts.get("min-support") {
        p.min_support_markers = v.parse().map_err(|_| "bad --min-support")?;
    }
    if let Some(v) = opts.get("min-coverage") {
        p.min_coverage = v.parse().map_err(|_| "bad --min-coverage (want 0..1)")?;
    }
    if let Some(v) = opts.get("min-abundance") {
        p.min_rel_abundance = v.parse().map_err(|_| "bad --min-abundance (want 0..1)")?;
    }
    // `--fixed-gate` turns both adaptations off (backwards compatible); the two halves can
    // also be controlled independently, which is what calibrating a large panel needs.
    if opts.contains_key("fixed-gate") {
        p.adaptive_singleton = false;
        p.adaptive_floor = false;
    }
    if opts.contains_key("no-adaptive-singleton") {
        p.adaptive_singleton = false;
    }
    if opts.contains_key("no-adaptive-floor") {
        p.adaptive_floor = false;
    }
    Ok(p)
}

fn cmd_profile(opts: &HashMap<String, String>) -> Result<(), String> {
    let db = StrainDb::load(Path::new(req(opts, "db")?)).map_err(|e| e.to_string())?;

    // Enzyme set: prefer the DB's recorded set (guarantees a match); else require --enzyme.
    let set: Vec<&Enzyme> = if !db.enzymes.is_empty() {
        parse_enzyme_set(&db.enzymes.join(",")).ok_or("DB records an unknown enzyme")?
    } else {
        enzyme_set(opts)?
    };

    let reads = PathBuf::from(req(opts, "reads")?);
    let counts = sample_marker_counts_stream(&reads, &set).map_err(|e| e.to_string())?;
    println!(
        "sample: {} distinct tag markers (enzymes: {}, threads: {})",
        counts.len(),
        enzyme_names(&set).join("+"),
        num_threads()
    );

    let params = parse_params(opts)?;
    let calls = profile(&db, &counts, &params);

    if calls.is_empty() {
        println!(
            "  (no strain/cluster resolved — insufficient strain-specific 2b tags for this \
             enzyme set; the species may still be present at Layer-1)"
        );
    } else {
        report(&calls);
    }

    if let Some(out) = opts.get("out") {
        write_pred_tsv(Path::new(out), &calls).map_err(|e| e.to_string())?;
        println!("predictions -> {out}");
    }
    Ok(())
}

/// Layer-1 outcome for one species in a multi-species sample.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SpeciesTier {
    /// Enough species-specific marker evidence to attempt within-species strain resolution.
    Resolved,
    /// Present, but below the evidence needed to resolve strains at this depth.
    DetectedNotResolved,
    /// Not enough evidence to call the species present.
    Absent,
}

/// Classify a species from ABSOLUTE species-specific marker evidence (never relative abundance).
///
/// `present` = its species-specific markers observed in the sample; `total` = its
/// species-specific markers that exist in the DB; `detect`/`floor` are absolute thresholds and
/// `frac` is the breadth fraction.
///
/// `reachable` is the fraction of the panel that *can* be observed at this species' estimated
/// depth (`1 − e^(−λ)`), and both gates are scaled by it. Without that scaling the fixed floor
/// of 200 markers is unreachable by construction in low-input or high-host samples: at 0.05×
/// depth only ~5% of any panel is visible, so a genuinely present species is filed as absent no
/// matter how clean the data is. The gate never falls below `detect`, which keeps random
/// sequencing errors — which produce scattered tags, not repeated hits on one species' panel —
/// from manufacturing calls. Pure + tested.
/// The adaptive relaxation is bounded: the floor is never scaled below this fraction of its
/// configured value.
///
/// Without a bound the scaling cancels out and stops being a gate at all. The observed count is
/// itself proportional to the reachable fraction (`present ≈ total·r` under Poisson sampling),
/// so testing `present ≥ floor·r` reduces to `total ≥ floor` — true for every species with a
/// panel above the floor, at every depth, leaving `detect` (10 markers) as the only real
/// threshold. Clamping the relaxation keeps a genuine floor at low depth while still recovering
/// species that an unscaled 200-marker bar makes unreachable by construction.
const MIN_FLOOR_FRACTION: f64 = 0.25;

fn species_tier(
    present: usize,
    total: usize,
    detect: usize,
    floor: usize,
    frac: f64,
    reachable: f64,
) -> SpeciesTier {
    let r = reachable.clamp(0.0, 1.0).max(MIN_FLOOR_FRACTION);
    let frac_gate = (frac.max(0.0) * total as f64 * r).ceil() as usize;
    let scaled_floor = (floor as f64 * r).ceil() as usize;
    let resolve_gate = scaled_floor.max(frac_gate).max(detect).max(1);
    let detect_gate = detect.min(resolve_gate);
    if present >= resolve_gate {
        SpeciesTier::Resolved
    } else if present >= detect_gate {
        SpeciesTier::DetectedNotResolved
    } else {
        SpeciesTier::Absent
    }
}

/// Per-species Layer-1 result carried out of the parallel map.
struct SpeciesResult {
    species: String,
    present_specific: usize,
    total_specific: usize,
    /// Estimated per-tag depth over this species' species-specific markers (zero-inclusive).
    lambda: f64,
    tier: SpeciesTier,
    calls: Vec<StrainCall>,
}

/// Multi-species strain profiling: digest the sample reads ONCE, then match the shared tag
/// counts against every per-species cluster DB in `--dbs <dir>`, in parallel across species.
/// This is the scalability advantage over running a full k-mer profiler once per species
/// (which re-counts k-mers every time).
fn cmd_multi_profile(opts: &HashMap<String, String>) -> Result<(), String> {
    let set = enzyme_set(opts)?;
    let dbs_dir = PathBuf::from(req(opts, "dbs")?);
    let reads = PathBuf::from(req(opts, "reads")?);

    // 1) digest sample reads ONCE (streamed: peak memory is one batch, not the whole file)
    let counts = sample_marker_counts_stream(&reads, &set).map_err(|e| e.to_string())?;

    // 2) collect + load per-species DBs
    let mut db_paths: Vec<PathBuf> = std::fs::read_dir(&dbs_dir)
        .map_err(|e| e.to_string())?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension().and_then(|x| x.to_str()) == Some("tsv")
                && p.file_name()
                    .and_then(|x| x.to_str())
                    .is_some_and(|n| !n.contains(".members."))
        })
        .collect();
    db_paths.sort();
    if db_paths.is_empty() {
        return Err("no *.tsv species DBs found in --dbs dir".into());
    }
    let mut loaded: Vec<(String, StrainDb)> = par_map(&db_paths, |path| {
        let sp = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        (sp, StrainDb::load(path).unwrap_or_default())
    });

    // 3) Layer-1 species gate (breadth-aware, three-tier). Strain markers are unique only
    //    *within* a species, so an absent species can be spuriously hit by a present relative's
    //    shared tags. We derive species-specific markers — tags carried by exactly ONE species
    //    across the panel (Strain2bScan's own species layer, same tag space as Fast2bRAD-M) —
    //    and decide per species from ABSOLUTE marker evidence, never relative abundance:
    //      total_specific   = the species' species-specific markers that exist in the DB
    //      present_specific = those observed in this sample (count >= 2)
    //    resolve_gate = max(min_species_markers, ceil(frac * total_specific)); the fraction term
    //    makes the bar comparable across species of very different panel sizes and maps to a
    //    minimum marker-panel coverage. detect_gate = min(min_species_detect, resolve_gate).
    //    Three outcomes: >= resolve_gate -> strain-resolved (Layer-2); >= detect_gate ->
    //    detected-not-resolvable (species-level only); else absent. All inputs come from
    //    Strain2bScan's own DBs + one digest of the raw reads — no external abundance needed.
    //    Defaults tuned on a 40-species panel (precision ~0.94 @ recall 1.0); calibrate per
    //    enzyme set and depth.
    let min_species_markers: usize = opts
        .get("min-species-markers")
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let min_species_marker_frac: f64 = opts
        .get("min-species-marker-frac")
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.0);
    let min_species_detect: usize = opts
        .get("min-species-detect")
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let mut species_degree: FxHashMap<Marker, u32> = FxHashMap::default();
    for (_, db) in &loaded {
        for &m in db.marker_degree.keys() {
            *species_degree.entry(m).or_insert(0) += 1;
        }
    }

    // Per species: the markers specific to it across the WHOLE panel (species degree == 1).
    let specific_sets: Vec<FxHashSet<Marker>> = loaded
        .iter()
        .map(|(_, db)| {
            db.marker_degree
                .keys()
                .copied()
                .filter(|m| species_degree.get(m).copied() == Some(1))
                .collect()
        })
        .collect();

    // Restrict detection AND quantification to those markers.
    //
    // Cluster-uniqueness is only defined within one species DB, so a tag can be unique to a
    // cluster here and still occur in a congener's genomes. When that congener is co-present —
    // S. aureus/S. epidermidis, the three streptococci, the two lactobacilli in MSA-1002, and
    // routinely in saliva — its reads land on that tag and inflate this cluster's depth. Using
    // the cross-species evidence for the species gate but not for quantification leaves exactly
    // that error in the abundances. `--no-cross-species-filter` disables this for comparison.
    let cross_species_filter = !opts.contains_key("no-cross-species-filter");
    if cross_species_filter {
        let (mut before, mut after) = (0usize, 0usize);
        for ((_, db), specific) in loaded.iter_mut().zip(&specific_sets) {
            before += db.marker_degree.len();
            db.restrict_to(specific);
            after += specific.len();
        }
        println!(
            "cross-species filter: {after}/{before} markers usable for quantification ({:.1}% shared with another species in the panel, excluded)",
            if before > 0 {
                100.0 * (before - after) as f64 / before as f64
            } else {
                0.0
            }
        );
    }
    let loaded = loaded;

    println!(
        "sample: {} distinct tag markers; {} species DBs; resolve-gate≥max({}, {:.0}%×panel), detect-gate≥{} (threads: {})",
        counts.len(),
        loaded.len(),
        min_species_markers,
        min_species_marker_frac * 100.0,
        min_species_detect,
        num_threads()
    );

    // 4) gate + strain-profile each species, in parallel.
    //    `--min-abundance` applies WITHIN a species, matching the primary output column: it
    //    decides which clusters of a present species are real. Whether the *species* is present
    //    at all is the Layer-1 species gate's decision (step 3), not this filter's.
    let params = parse_params(opts)?;
    let order: Vec<usize> = (0..loaded.len()).collect();
    let per_species: Vec<SpeciesResult> = par_map(&order, |&i| {
        let (species, db) = &loaded[i];
        let specific = &specific_sets[i];
        let total_specific = specific.len();
        // Zero-inclusive mean over the species-specific panel — an absolute depth estimate
        // that does not presuppose the species passed any gate (so it is not circular).
        let observed: u64 = specific
            .iter()
            .map(|m| counts.get(m).copied().unwrap_or(0) as u64)
            .sum();
        let lambda = if total_specific > 0 {
            observed as f64 / total_specific as f64
        } else {
            0.0
        };
        let min_count = if params.adaptive_singleton {
            strain2bscan::identify::min_count_for(lambda)
        } else {
            2
        };
        let present_specific = specific
            .iter()
            .filter(|m| counts.get(*m).copied().unwrap_or(0) >= min_count)
            .count();
        let reachable = if params.adaptive_floor {
            detectable_fraction(lambda)
        } else {
            1.0
        };
        let tier = species_tier(
            present_specific,
            total_specific,
            min_species_detect,
            min_species_markers,
            min_species_marker_frac,
            reachable,
        );
        let calls = if tier == SpeciesTier::Resolved {
            profile(db, &counts, &params)
        } else {
            Vec::new()
        };
        SpeciesResult {
            species: species.clone(),
            present_specific,
            total_specific,
            lambda,
            tier,
            calls,
        }
    });

    // 5) Quantification. Two scopes, both reported:
    //
    //    `abundance`        — WITHIN-species (sums to 1.0 per species). This is the primary
    //                         number and the one Layer-2 is actually answering: given that this
    //                         species is present, how is it split across strains/clusters?
    //                         Species-level abundance is Fast2bRAD-M's job, not this layer's.
    //    `global_abundance` — cross-species, derived from the absolute `depth`. Provided so a
    //                         community composition can be plotted without re-deriving it:
    //                         per-species fractions cannot be concatenated into one, because
    //                         each species sums to 1.0 independently and a 10x-more-abundant
    //                         species would otherwise look identical to a rare one.
    //
    //    `depth` is reads-per-tag on single-copy markers (one copy per cell), so it is
    //    comparable between species and these are cell fractions.
    //
    //    IMPORTANT: the denominator spans only the clusters actually emitted. Species that were
    //    detected but not strain-resolvable, and clusters dropped by --min-abundance, carry real
    //    biomass that is NOT in it, so `global_abundance` is a composition of the strain-resolved
    //    fraction of the sample, not of the whole community. The summary line below reports how
    //    many species fell outside it; use `depth` directly if you need absolute quantities.
    let depth_sum: f64 = per_species
        .iter()
        .flat_map(|r| r.calls.iter())
        .map(|c| c.depth)
        .sum();
    let global_of = |c: &StrainCall| {
        if depth_sum > 0.0 {
            c.depth / depth_sum
        } else {
            0.0
        }
    };

    // `sample_fraction` — share of ALL tag observations in the sample.
    //
    // `global_abundance`'s denominator is whatever this run happened to resolve, so it is not
    // comparable between samples: change the depth, the reference panel, or the number of
    // unknown organisms present, and the denominator moves. In a defined mock (membership known,
    // everything detectable) that is harmless; in a real sample, where most of the community may
    // have no reference at all, it silently rescales everything.
    //
    // This column uses a sample-intrinsic denominator instead — the total tag observations —
    // which is fixed by the sequencing, not by how well profiling went. A cluster's numerator is
    // its estimated tag mass, `depth * n_markers` (per-tag depth times how many tags it carries;
    // the units of `depth` alone are reads-per-tag, which cannot be compared against a read
    // total). Clusters sharing a tag each contribute their own depth to it, so the shares add up
    // without double counting. The remainder is reported as unclassified.
    let total_tags: u64 = counts.values().map(|&c| c as u64).sum();
    let mass_of = |c: &StrainCall| c.depth * c.n_markers as f64;
    let sample_fraction_of = |c: &StrainCall| {
        if total_tags > 0 {
            mass_of(c) / total_tags as f64
        } else {
            0.0
        }
    };

    // 6) report, grouped by species, most abundant cluster first within each
    let mut flat: Vec<(&str, &StrainCall, f64, f64)> = per_species
        .iter()
        .flat_map(|r| {
            r.calls
                .iter()
                .map(move |c| (r.species.as_str(), c, global_of(c), sample_fraction_of(c)))
        })
        .collect();
    flat.sort_by(|a, b| {
        a.0.cmp(b.0).then(
            b.1.rel_abundance
                .partial_cmp(&a.1.rel_abundance)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
    });

    // Append-only column order: the first six match what `multi-profile` wrote before the
    // accuracy rewrite, so positional readers keep working; new columns go on the end.
    println!(
        "#species\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction"
    );
    for (species, c, g, sf) in &flat {
        println!(
            "  {}\t{}\t{:.6}\t{:.2}\t{:.0}\t{:.3}\t{:.6}\t{:.6}",
            species, c.name, c.rel_abundance, c.coverage, c.support, c.depth, g, sf
        );
    }

    let (mut n_resolved, mut n_detected) = (0usize, 0usize);
    for r in &per_species {
        match r.tier {
            SpeciesTier::Resolved => {
                n_resolved += 1;
                if r.calls.is_empty() {
                    println!(
                        "  {}\t[strain-resolved, no cluster above threshold]\tmarkers={}/{} (depth {:.2}x)",
                        r.species, r.present_specific, r.total_specific, r.lambda
                    );
                }
            }
            SpeciesTier::DetectedNotResolved => {
                n_detected += 1;
                let breadth = if r.total_specific > 0 {
                    100.0 * r.present_specific as f64 / r.total_specific as f64
                } else {
                    0.0
                };
                println!(
                    "  {}\t[detected, not strain-resolvable]\tmarkers={}/{} ({:.1}%, depth {:.2}x)",
                    r.species, r.present_specific, r.total_specific, breadth, r.lambda
                );
            }
            SpeciesTier::Absent => {}
        }
    }
    println!(
        "summary: {}/{} species strain-resolved ({} strain calls), {} detected-not-resolvable, {} absent",
        n_resolved,
        loaded.len(),
        flat.len(),
        n_detected,
        loaded.len() - n_resolved - n_detected
    );
    // How much of the sequencing the strain calls actually account for. The remainder is
    // unresolved species, organisms with no reference in the panel, host DNA and sequencing
    // error — real signal that `global_abundance`'s denominator silently omits.
    let classified: f64 = flat.iter().map(|(_, c, _, _)| mass_of(c)).sum();
    if total_tags > 0 {
        let pct = 100.0 * classified / total_tags as f64;
        println!(
            "coverage of sample: strain calls account for {:.1}% of {} tag observations ({:.1}% unclassified — unresolved species, no reference, host, error)",
            pct.min(100.0),
            total_tags,
            (100.0 - pct).max(0.0)
        );
    }

    if let Some(out) = opts.get("out") {
        use std::io::Write;
        let mut w = std::fs::File::create(out).map_err(|e| e.to_string())?;
        writeln!(
            w,
            "#species\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction\tn_markers"
        )
        .map_err(|e| e.to_string())?;
        for (species, c, g, sf) in &flat {
            writeln!(
                w,
                "{}\t{}\t{:.6}\t{:.4}\t{:.0}\t{:.4}\t{:.6}\t{:.6}\t{}",
                species, c.name, c.rel_abundance, c.coverage, c.support, c.depth, g, sf, c.n_markers
            )
            .map_err(|e| e.to_string())?;
        }
        println!("predictions -> {out}");
    }
    Ok(())
}

fn cmd_info(opts: &HashMap<String, String>) -> Result<(), String> {
    let db = StrainDb::load(Path::new(req(opts, "db")?)).map_err(|e| e.to_string())?;
    println!(
        "enzymes: {}",
        if db.enzymes.is_empty() {
            "(unspecified)".into()
        } else {
            db.enzymes.join("+")
        }
    );
    print_stats(&db);
    for (i, name) in db.strain_names.iter().enumerate() {
        println!(
            "  [{i}] {name}: {} markers ({} unique)",
            db.strain_markers[i].len(),
            db.unique_marker_count(i)
        );
    }
    Ok(())
}

/// In-memory conspecific demo (no files needed): shows the ported Layer-2 resolving a
/// 70/30 mixture of 2 near-identical strains while the naive scorer over-calls all 4.
fn cmd_demo() -> Result<(), String> {
    let core: Vec<Marker> = (0..200).collect();
    let mut strains = Vec::new();
    for s in 0..4u64 {
        let mut m = core.clone();
        m.extend((0..50).map(|i| 1_000_000 + s * 50 + i));
        strains.push((format!("strain{s}"), m));
    }
    let db = StrainDb::build(strains);

    let mixture = [(0usize, 0.7f64), (2, 0.3)];
    let mut present: HashMap<Marker, f64> = HashMap::new();
    for &(j, ab) in &mixture {
        for &m in &db.strain_markers[j] {
            *present.entry(m).or_insert(0.0) += ab;
        }
    }
    let mut counts = MarkerCounts::default();
    for (m, frac) in present {
        let c = (30.0 * frac).round() as u32;
        if c > 0 {
            counts.insert(m, c);
        }
    }
    for e in 0..50u64 {
        counts.insert(9_000_000 + e, 1);
    }

    println!("== Demo: 4 conspecific strains (200 shared + 50 private each) ==");
    println!("truth: strain0=0.70, strain2=0.30\n");
    println!("[ported StrainScan Layer-2]");
    report(&profile(&db, &counts, &Params::default()));
    let naive = naive_profile(&db, &counts, 1240.0);
    println!(
        "\n[naive strainscan-rust-style scoring]  -> calls {} strains: {:?}  (over-call: shared core alone clears the threshold)",
        naive.len(),
        naive.iter().map(|&j| db.strain_names[j].clone()).collect::<Vec<_>>()
    );
    Ok(())
}

/// CST demo: a species with 4 genomes forming 2 clusters. Shows single-linkage clustering,
/// marker classification, and cluster-resolution profiling of a 70/30 cluster mixture.
fn cmd_cst_demo() -> Result<(), String> {
    let core: Vec<Marker> = (0..200).collect();
    let clu_a: Vec<Marker> = (200..240).collect();
    let clu_b: Vec<Marker> = (300..340).collect();
    let mk = |name: &str, extra: &[Marker], base: Marker| {
        let mut v = core.clone();
        v.extend_from_slice(extra);
        v.extend((0..3).map(|i| base + i));
        (name.to_string(), v.clone(), v) // (name, single-copy, full) — no multi-copy here
    };
    let genomes = vec![
        mk("g0", &clu_a, 1000),
        mk("g1", &clu_a, 1100),
        mk("g2", &clu_b, 2000),
        mk("g3", &clu_b, 2100),
    ];

    println!("== CST demo: 1 species, 4 genomes (g0/g1 ~identical, g2/g3 ~identical) ==");
    let cst = SpeciesCst::build(genomes, DEFAULT_SIMILARITY, false);
    println!("single-linkage @ 0.95 -> {} clusters:", cst.n_clusters());
    for (cid, members) in cst.clusters.iter().enumerate() {
        let names: Vec<&str> = members
            .iter()
            .map(|&g| cst.genome_names[g].as_str())
            .collect();
        println!("  C{cid}: {}", names.join(", "));
    }
    let s = cst.marker_class_summary();
    println!(
        "marker classes: species_core={} cluster_specific={} strain_specific={}",
        s.get("species_core").unwrap_or(&0),
        s.get("cluster_specific").unwrap_or(&0),
        s.get("strain_specific").unwrap_or(&0),
    );

    let db = cst.cluster_db();
    let mut present: HashMap<Marker, f64> = HashMap::new();
    for (cid, ab) in [(0usize, 0.7f64), (1, 0.3)] {
        for &m in &db.strain_markers[cid] {
            *present.entry(m).or_insert(0.0) += ab;
        }
    }
    let mut counts = MarkerCounts::default();
    for (m, frac) in present {
        let c = (30.0 * frac).round() as u32;
        if c > 0 {
            counts.insert(m, c);
        }
    }
    println!("\nprofiling cluster mixture truth C0=0.70, C1=0.30:");
    report(&profile(&db, &counts, &Params::default()));
    Ok(())
}

fn report(calls: &[StrainCall]) {
    if calls.is_empty() {
        println!("  (no strains passed thresholds)");
        return;
    }
    for c in calls {
        println!(
            "  {:<12} abundance={:>6.2}%  depth={:>7.3}x  coverage={:>6.2}%  support={:.0}",
            c.name,
            c.rel_abundance * 100.0,
            c.depth,
            c.coverage * 100.0,
            c.support
        );
    }
}

/// Write predictions as `name<TAB>abundance<TAB>coverage<TAB>support<TAB>depth<TAB>n_markers`.
///
/// Column order is **append-only**: the first five columns are exactly those written before the
/// accuracy rewrite, so scripts that index positionally keep working. Inserting `depth` in the
/// middle silently shifts `coverage` and `support` by one, which does not fail loudly — it just
/// makes every downstream number wrong.
///
/// `depth` is the absolute per-tag read depth; unlike `abundance` (normalized within this DB)
/// it is comparable across separately-profiled species. `n_markers` is the cluster's tag count,
/// needed to convert depth into a DNA-mass share.
fn write_pred_tsv(path: &Path, calls: &[StrainCall]) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::fs::File::create(path)?;
    writeln!(w, "#cluster\tabundance\tcoverage\tsupport\tdepth\tn_markers")?;
    for c in calls {
        writeln!(
            w,
            "{}\t{:.6}\t{:.4}\t{:.0}\t{:.4}\t{}",
            c.name, c.rel_abundance, c.coverage, c.support, c.depth, c.n_markers
        )?;
    }
    Ok(())
}

fn print_stats(db: &StrainDb) {
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
    use super::{species_tier, SpeciesTier};
    use strain2bscan::identify::detectable_fraction;

    /// Ample depth: essentially the whole panel is reachable, so the gates are the original
    /// absolute ones.
    const FULL: f64 = 1.0;

    #[test]
    fn absolute_floor_gates_when_no_fraction() {
        // frac = 0 -> resolve gate is the absolute floor (200); detect gate = min(10, 200) = 10.
        assert_eq!(species_tier(250, 5000, 10, 200, 0.0, FULL), SpeciesTier::Resolved);
        assert_eq!(species_tier(50, 5000, 10, 200, 0.0, FULL), SpeciesTier::DetectedNotResolved);
        assert_eq!(species_tier(5, 5000, 10, 200, 0.0, FULL), SpeciesTier::Absent);
    }

    #[test]
    fn breadth_fraction_raises_the_bar_for_large_panels() {
        // 10% of a 5000-marker panel = 500 > floor 200, so 300 observed is below the resolve gate
        // even though it clears the absolute floor. This is the whole point of the breadth term.
        assert_eq!(species_tier(300, 5000, 10, 200, 0.10, FULL), SpeciesTier::DetectedNotResolved);
        assert_eq!(species_tier(600, 5000, 10, 200, 0.10, FULL), SpeciesTier::Resolved);
    }

    #[test]
    fn small_panel_species_can_still_be_detected() {
        // total 150 < floor 200 -> can never clear the resolve gate, but the detect gate
        // (min(10, 200) = 10) still flags presence rather than dropping it silently.
        assert_eq!(species_tier(150, 150, 10, 200, 0.0, FULL), SpeciesTier::DetectedNotResolved);
        assert_eq!(species_tier(5, 150, 10, 200, 0.0, FULL), SpeciesTier::Absent);
    }

    /// At low depth only a few percent of any panel is observable, so an unscaled 200-marker
    /// floor is unreachable by construction and a genuinely present species is filed as
    /// unresolvable. The floor relaxes — but only to `MIN_FLOOR_FRACTION` of its configured
    /// value (200 -> 50), never all the way down.
    #[test]
    fn low_depth_relaxes_the_floor_but_only_to_the_bound() {
        let reachable = detectable_fraction(0.05);
        assert!(reachable < 0.05);
        // 60 observed markers: below the full-depth floor of 200 ...
        assert_eq!(species_tier(60, 5000, 10, 200, 0.0, FULL), SpeciesTier::DetectedNotResolved);
        // ... but above the relaxed floor of 50 at low depth.
        assert_eq!(species_tier(60, 5000, 10, 200, 0.0, reachable), SpeciesTier::Resolved);
    }

    /// The relaxation must be **bounded**. Scaling the floor by the reachable fraction `r` alone
    /// is self-cancelling: the observed count is itself ~ `total * r`, so `present >= floor * r`
    /// reduces to `total >= floor` — independent of depth. That leaves `detect` (10 markers) as
    /// the only real threshold, and a 20 000-marker panel hit by 40 stray tags would be declared
    /// strain-resolvable. The bound keeps a genuine floor.
    #[test]
    fn adaptive_floor_does_not_cancel_itself_away() {
        let r = detectable_fraction(0.002); // ~0.2% reachable
        // 40 stray tags on a 20k panel: detected as present, but NOT strain-resolvable.
        assert_eq!(
            species_tier(40, 20_000, 10, 200, 0.0, r),
            SpeciesTier::DetectedNotResolved
        );
        // The bounded floor is 50, so real evidence still resolves.
        assert_eq!(species_tier(50, 20_000, 10, 200, 0.0, r), SpeciesTier::Resolved);
    }

    /// The gate must never fall below the absolute detect floor, at any depth.
    #[test]
    fn adaptive_gate_never_falls_below_the_detect_floor() {
        for lambda in [0.0, 1e-6, 0.001, 0.01] {
            let r = detectable_fraction(lambda);
            assert_eq!(species_tier(9, 5000, 10, 200, 0.0, r), SpeciesTier::Absent);
            assert_eq!(
                species_tier(10, 5000, 10, 200, 0.0, r),
                SpeciesTier::DetectedNotResolved
            );
        }
    }
}
