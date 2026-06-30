//! strain2bscan CLI (prototype).
//!
//!   strain2bscan build    --genomes <dir> --enzyme <set> --out <db.tsv> [--max-contigs N] [--min-tag-fraction F]
//!   strain2bscan cluster  --genomes <dir> --enzyme <set> --out <clusterdb.tsv> [--similarity 0.95] [--max-contigs N] [--min-tag-fraction F]
//!   strain2bscan profile  --db <db.tsv> --reads <fastx> [--enzyme <set>] [--out pred.tsv] [--min-support N]
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
use strain2bscan::identify::{naive_profile, profile, Params, StrainCall};
use strain2bscan::markers::{
    genome_marker_counts_multi, read_fastx, sample_marker_counts_multi_par, single_copy_markers,
    Marker,
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
                 strain2bscan cluster  --genomes <dir> --enzyme <set> --out <clusterdb.tsv> [--similarity 0.95] [--max-contigs N] [--min-tag-fraction F]\n  \
                 strain2bscan profile  --db <db.tsv> --reads <fastx> [--enzyme <set>] [--out pred.tsv] [--min-support N]\n  \
                 strain2bscan multi-profile --dbs <dir> --reads <fastx> --enzyme <set>   (many species, sample digested once)\n  \
                 strain2bscan info     --db <db.tsv>\n  \
                 strain2bscan evaluate --pred <pred.tsv> --truth <truth.tsv> [--present 0.01]\n  \
                 strain2bscan demo | cst-demo\n\n\
                 <set> = all | BcgI | BcgI,CspCI  (use BcgI for 2bRAD data; all for conventional metagenomes)"
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
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
        if matches!(ext, "fa" | "fasta" | "fna") {
            paths.push(path);
        }
    }
    if paths.is_empty() {
        return Err("no FASTA genomes (.fa/.fasta/.fna) found".into());
    }
    paths.sort(); // deterministic genome order regardless of threading
    let results: Vec<Result<GenomeRec, String>> = par_map(&paths, |path| {
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        let seqs = read_fastx(path).map_err(|e| e.to_string())?;
        let n_contigs = seqs.len();
        let counts = genome_marker_counts_multi(&seqs, enzymes);
        Ok(GenomeRec { name, n_contigs, markers: single_copy_markers(&counts) })
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
) -> Result<Vec<(String, Vec<Marker>)>, String> {
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
    Ok(rep.kept.into_iter().map(|g| (g.name, g.markers)).collect())
}

fn cmd_build(opts: &HashMap<String, String>) -> Result<(), String> {
    let set = enzyme_set(opts)?;
    let genomes = PathBuf::from(req(opts, "genomes")?);
    let out = PathBuf::from(req(opts, "out")?);

    let strains = digest_and_filter(&genomes, &set, opts)?;
    for (name, m) in &strains {
        println!("  {name}: {} single-copy tag markers", m.len());
    }
    let mut db = StrainDb::build(strains);
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

    let strains = digest_and_filter(&genomes, &set, opts)?;
    let n_genomes = strains.len();
    let cst = SpeciesCst::build(strains, similarity);
    let method = if n_genomes > MINHASH_ABOVE {
        "MinHash"
    } else {
        "exact-Jaccard"
    };
    println!(
        "clustered {} genomes into {} cluster(s) @ similarity {similarity} (enzymes: {}, threads: {}, clustering: {})",
        cst.genome_names.len(),
        cst.n_clusters(),
        enzyme_names(&set).join("+"),
        num_threads(),
        method
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

fn cmd_profile(opts: &HashMap<String, String>) -> Result<(), String> {
    let db = StrainDb::load(Path::new(req(opts, "db")?)).map_err(|e| e.to_string())?;

    // Enzyme set: prefer the DB's recorded set (guarantees a match); else require --enzyme.
    let set: Vec<&Enzyme> = if !db.enzymes.is_empty() {
        parse_enzyme_set(&db.enzymes.join(",")).ok_or("DB records an unknown enzyme")?
    } else {
        enzyme_set(opts)?
    };

    let reads = PathBuf::from(req(opts, "reads")?);
    let seqs = read_fastx(&reads).map_err(|e| e.to_string())?;
    let counts = sample_marker_counts_multi_par(&seqs, &set);
    println!(
        "sample: {} distinct tag markers (enzymes: {}, threads: {})",
        counts.len(),
        enzyme_names(&set).join("+"),
        num_threads()
    );

    let mut params = Params::default();
    if let Some(v) = opts.get("min-support") {
        params.min_support_markers = v.parse().map_err(|_| "bad --min-support")?;
    }
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

/// Multi-species strain profiling: digest the sample reads ONCE, then match the shared tag
/// counts against every per-species cluster DB in `--dbs <dir>`, in parallel across species.
/// This is the scalability advantage over running a full k-mer profiler once per species
/// (which re-counts k-mers every time).
fn cmd_multi_profile(opts: &HashMap<String, String>) -> Result<(), String> {
    let set = enzyme_set(opts)?;
    let dbs_dir = PathBuf::from(req(opts, "dbs")?);
    let reads = PathBuf::from(req(opts, "reads")?);

    // 1) digest sample reads ONCE
    let seqs = read_fastx(&reads).map_err(|e| e.to_string())?;
    let counts = sample_marker_counts_multi_par(&seqs, &set);

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
    let loaded: Vec<(String, StrainDb)> = par_map(&db_paths, |path| {
        let sp = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        (sp, StrainDb::load(path).unwrap_or_default())
    });

    // 3) Layer-1 species gate. Strain markers are unique only *within* a species, so an
    //    absent species can be spuriously hit by a present relative's shared tags. Derive
    //    species-specific markers (carried by exactly ONE species across the panel — the
    //    Fast2bRAD species layer) and require enough of them present before profiling strains.
    // Default tuned on a 40-species real-genome panel (precision ~0.94 @ recall 1.0); the
    // proper production gate is Fast2bRAD-M's species-level call. Calibrate per panel/depth.
    let min_species_markers: usize = opts
        .get("min-species-markers")
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let mut species_degree: HashMap<Marker, u32> = HashMap::new();
    for (_, db) in &loaded {
        for &m in db.marker_degree.keys() {
            *species_degree.entry(m).or_insert(0) += 1;
        }
    }

    println!(
        "sample: {} distinct tag markers; {} species DBs; species-gate≥{} specific markers (threads: {})",
        counts.len(),
        loaded.len(),
        min_species_markers,
        num_threads()
    );

    // 4) gate + strain-profile each species, in parallel
    let params = Params::default();
    let per_species: Vec<(String, usize, Vec<StrainCall>)> = par_map(&loaded, |(species, db)| {
        let present_specific = db
            .marker_degree
            .keys()
            .filter(|m| {
                species_degree.get(m).copied().unwrap_or(0) == 1
                    && counts.get(m).copied().unwrap_or(0) >= 2
            })
            .count();
        if present_specific < min_species_markers {
            return (species.clone(), present_specific, Vec::new());
        }
        (
            species.clone(),
            present_specific,
            profile(db, &counts, &params),
        )
    });

    let mut total_calls = 0;
    for (species, _, calls) in &per_species {
        total_calls += calls.len();
        for c in calls {
            println!(
                "  {species}\t{}\t{:.4}\t{:.2}\t{:.0}",
                c.name, c.rel_abundance, c.coverage, c.support
            );
        }
    }
    println!(
        "detected strains in {}/{} species ({} strain calls total)",
        per_species.iter().filter(|(_, _, c)| !c.is_empty()).count(),
        loaded.len(),
        total_calls
    );
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
    let mut counts: HashMap<Marker, u32> = HashMap::new();
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
    let mk = |extra: &[Marker], base: Marker| {
        let mut v = core.clone();
        v.extend_from_slice(extra);
        v.extend((0..3).map(|i| base + i));
        v
    };
    let genomes = vec![
        ("g0".to_string(), mk(&clu_a, 1000)),
        ("g1".to_string(), mk(&clu_a, 1100)),
        ("g2".to_string(), mk(&clu_b, 2000)),
        ("g3".to_string(), mk(&clu_b, 2100)),
    ];

    println!("== CST demo: 1 species, 4 genomes (g0/g1 ~identical, g2/g3 ~identical) ==");
    let cst = SpeciesCst::build(genomes, DEFAULT_SIMILARITY);
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
    let mut counts: HashMap<Marker, u32> = HashMap::new();
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
            "  {:<12} abundance={:>6.2}%  coverage={:>6.2}%  support={:.0}",
            c.name,
            c.rel_abundance * 100.0,
            c.coverage * 100.0,
            c.support
        );
    }
}

/// Write predictions as `name<TAB>abundance<TAB>coverage<TAB>support` (header commented).
fn write_pred_tsv(path: &Path, calls: &[StrainCall]) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::fs::File::create(path)?;
    writeln!(w, "#cluster\tabundance\tcoverage\tsupport")?;
    for c in calls {
        writeln!(
            w,
            "{}\t{:.6}\t{:.4}\t{:.0}",
            c.name, c.rel_abundance, c.coverage, c.support
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
