//! Build a strain DB from genomes.

use crate::cli::BuildArgs;
use crate::{digest_and_filter, kmer_params, marker_source_for_build, print_stats};
use strain2bscan::db::StrainDb;

pub fn run(args: &BuildArgs) -> Result<(), String> {
    let (k, scale) = kmer_params(Some(args.kmer_size), Some(args.sketch_scale))?;
    let source = marker_source_for_build(&args.marker_source, args.enzyme.as_deref(), k, scale)?;
    let genomes = &args.genomes;
    let out = &args.out;

    let recs = digest_and_filter(genomes, &source, args.max_contigs, args.min_tag_fraction)?;
    for r in &recs {
        let what = match source {
            crate::MarkerSource::Enzyme(_) => "tag markers",
            crate::MarkerSource::Kmer { .. } => "k-mer markers",
        };
        println!("  {}: {} single-copy {}", r.name, r.markers.len(), what);
    }
    let mut db = StrainDb::build(recs.into_iter().map(|r| (r.name, r.markers)).collect());
    db.enzymes = source.db_token();
    db.save(out).map_err(|e| e.to_string())?;
    print_stats(&db);
    println!("saved DB ({}) -> {}", db.enzymes.join("+"), out.display());
    Ok(())
}
