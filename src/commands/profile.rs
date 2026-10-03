//! Profile a single species sample against a strain DB.

use crate::cli::ProfileArgs;
use crate::commands::parse_params;
use crate::report::{reads_path, report, write_pred_tsv};
use crate::{resolve_sample_source, MarkerSource};
use strain2bscan::db::StrainDb;
use strain2bscan::identify::{resolve_layer1, tree_utility, Layer1};
use strain2bscan::parallel::num_threads;

pub fn run(args: &ProfileArgs) -> Result<(), String> {
    let db = StrainDb::load(&args.db).map_err(|e| e.to_string())?;
    let source = resolve_sample_source(
        &db,
        args.marker_source.as_deref(),
        args.enzyme.as_deref(),
        args.kmer_size,
        args.sketch_scale,
    )?;

    let reads = reads_path(args.reads.to_str().unwrap_or(""))?;
    let counts = source.sample_counts(&reads)?;
    println!(
        "sample: {} distinct {} ({}, threads: {})",
        counts.len(),
        match source {
            MarkerSource::Enzyme(_) => "tag markers",
            MarkerSource::Kmer { .. } => "k-mer markers",
        },
        source.describe(),
        num_threads()
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
        false,
        false,
    )?;

    let chosen = resolve_layer1(&db, &params);
    if params.layer1 == Layer1::Auto {
        match tree_utility(&db, params.min_support_markers) {
            Some(u) => println!(
                "layer1: {} (auto — {} cluster(s) below the support floor, {} informative internal node(s), {} fallback node(s))",
                if chosen == Layer1::Cst { "cst" } else { "unique" },
                u.rescuable,
                u.informative_nodes,
                u.fallback_nodes
            ),
            None => println!("layer1: unique (auto — this DB carries no tree)"),
        }
    }
    let calls = strain2bscan::identify::profile(&db, &counts, &params);

    if calls.is_empty() {
        match source {
            MarkerSource::Enzyme(_) => println!(
                "  (no strain/cluster resolved — insufficient strain-specific 2b tags for this \
                 enzyme set; the species may still be present at Layer-1)"
            ),
            MarkerSource::Kmer { .. } => println!(
                "  (no strain/cluster resolved — insufficient strain-specific k-mer markers for \
                 this sketch; the species may still be present at Layer-1)"
            ),
        }
    } else {
        report(&calls);
    }

    if let Some(out) = &args.out {
        write_pred_tsv(out, &calls).map_err(|e| e.to_string())?;
        println!("predictions -> {}", out.display());
    }
    Ok(())
}
