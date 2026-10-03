//! Batch multi-species profiling from a manifest.

use std::io::Write;

use crate::cli::BatchArgs;
use crate::commands::parse_params;
use crate::panel::{load_panel, profile_sample, species_gate};
use crate::report::read_manifest;

pub fn run(args: &BatchArgs) -> Result<(), String> {
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
    let samples = read_manifest(&args.manifest)?;

    let mut w =
        std::io::BufWriter::new(std::fs::File::create(&args.out).map_err(|e| e.to_string())?);
    writeln!(
        w,
        "#sample\tspecies\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction\tn_markers"
    )
    .map_err(|e| e.to_string())?;

    for (i, s) in samples.iter().enumerate() {
        eprintln!(
            "[batch {}/{}] sample {}: {}",
            i + 1,
            samples.len(),
            s.name,
            match &s.reads2 {
                Some(r2) => format!("{} + {}", s.reads1.display(), r2.display()),
                None => s.reads1.display().to_string(),
            }
        );
        let counts = panel
            .source
            .sample_counts_paired(&s.reads1, s.reads2.as_deref())
            .map_err(|e| format!("sample '{}': {e}", s.name))?;
        let result = profile_sample(&panel, &counts, &gate, &params, args.min_global_abundance)?;
        eprintln!(
            "[batch {}/{}] sample {}: {} distinct markers, {} strain call(s)",
            i + 1,
            samples.len(),
            s.name,
            counts.len(),
            result.rows.len()
        );
        for r in &result.rows {
            let c = &r.call;
            writeln!(
                w,
                "{}\t{}\t{}\t{:.6}\t{:.4}\t{:.0}\t{:.4}\t{:.6}\t{:.6}\t{}",
                s.name,
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
    }
    w.flush().map_err(|e| e.to_string())?;
    println!(
        "predictions ({} samples) -> {}",
        samples.len(),
        args.out.display()
    );
    Ok(())
}
