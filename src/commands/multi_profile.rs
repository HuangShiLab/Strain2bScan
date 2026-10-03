//! Multi-species strain profiling against per-species DBs.

use crate::cli::MultiProfileArgs;
use crate::commands::parse_params;
use crate::panel::{load_panel, profile_sample, species_gate, SpeciesTier};
use crate::report::reads_path;
use crate::MarkerSource;
use strain2bscan::identify::Layer1;
use strain2bscan::parallel::num_threads;

pub fn run(args: &MultiProfileArgs) -> Result<(), String> {
    let reads = reads_path(args.reads.to_str().unwrap_or(""))?;
    let panel = load_panel(
        &args.dbs,
        args.marker_source.as_deref(),
        args.enzyme.as_deref(),
        args.kmer_size,
        args.sketch_scale,
        args.no_cross_species_filter,
    )?;
    let gate = species_gate(
        args.min_species_markers,
        args.min_species_marker_frac,
        args.min_species_detect,
    );
    let params = parse_params(
        args.min_support,
        args.min_coverage,
        args.min_abundance,
        args.trace_gap,
        args.trace_floor,
        args.layer1.as_deref(),
        args.layer2.as_deref(),
        args.enet_alpha,
        args.min_consistency,
        args.fixed_gate,
        args.no_adaptive_singleton,
        args.no_adaptive_floor,
    )?;

    let counts = panel.source.sample_counts(&reads)?;

    match &panel.source {
        MarkerSource::Enzyme(_) => println!(
            "sample: {} distinct tag markers; {} species DBs; resolve-gate≥max({}, {:.0}%×panel), detect-gate≥{} (threads: {})",
            counts.len(),
            panel.loaded.len(),
            gate.min_markers,
            gate.min_frac * 100.0,
            gate.min_detect,
            num_threads()
        ),
        MarkerSource::Kmer { .. } => println!(
            "sample: {} distinct k-mer markers ({}); {} species DBs; resolve-gate≥max({}, {:.0}%×panel), detect-gate≥{} (threads: {})",
            counts.len(),
            panel.source.describe(),
            panel.loaded.len(),
            gate.min_markers,
            gate.min_frac * 100.0,
            gate.min_detect,
            num_threads()
        ),
    }

    let result = profile_sample(&panel, &counts, &gate, &params, args.min_global_abundance)?;

    println!(
        "#species\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction"
    );
    for r in &result.rows {
        let c = &r.call;
        println!(
            "  {}\t{}\t{:.6}\t{:.2}\t{:.0}\t{:.3}\t{:.6}\t{:.6}",
            r.species,
            c.name,
            c.rel_abundance,
            c.coverage,
            c.support,
            c.depth,
            r.global_abundance,
            r.sample_fraction
        );
    }

    let (mut n_resolved, mut n_detected) = (0usize, 0usize);
    for r in &result.per_species {
        match r.tier {
            SpeciesTier::Resolved => {
                n_resolved += 1;
                if params.layer1 == Layer1::Auto {
                    match r.tree {
                        Some(u) => println!(
                            "  {}\tlayer1={} (auto — {} cluster(s) below the support floor, {} informative internal node(s), {} fallback node(s))",
                            r.species,
                            if r.layer1 == Layer1::Cst { "cst" } else { "unique" },
                            u.rescuable,
                            u.informative_nodes,
                            u.fallback_nodes
                        ),
                        None => println!(
                            "  {}\tlayer1=unique (auto — this DB carries no tree)",
                            r.species
                        ),
                    }
                }
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
        panel.loaded.len(),
        result.rows.len(),
        n_detected,
        panel.loaded.len() - n_resolved - n_detected
    );

    let mass_of = |c: &strain2bscan::identify::StrainCall| c.depth * c.n_markers as f64;
    let classified: f64 = result.rows.iter().map(|r| mass_of(&r.call)).sum();
    if result.total_tags > 0 {
        let pct = 100.0 * classified / result.total_tags as f64;
        println!(
            "coverage of sample: strain calls account for {:.1}% of {} tag observations ({:.1}% unclassified — unresolved species, no reference, host, error)",
            pct.min(100.0),
            result.total_tags,
            (100.0 - pct).max(0.0)
        );
    }

    if let Some(out) = &args.out {
        use std::io::Write;
        let mut w = std::fs::File::create(out).map_err(|e| e.to_string())?;
        writeln!(
            w,
            "#species\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction\tn_markers"
        )
        .map_err(|e| e.to_string())?;
        for r in &result.rows {
            let c = &r.call;
            writeln!(
                w,
                "{}\t{}\t{:.6}\t{:.4}\t{:.0}\t{:.4}\t{:.6}\t{:.6}\t{}",
                r.species,
                c.name,
                c.rel_abundance,
                c.coverage,
                c.support,
                c.depth,
                r.global_abundance,
                r.sample_fraction,
                c.n_markers
            )
            .map_err(|e| e.to_string())?;
        }
        println!("predictions -> {}", out.display());
    }
    Ok(())
}
