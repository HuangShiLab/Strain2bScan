//! Measure whether a Cluster Search Tree could work on a panel.

use crate::cli::DiagnoseTreeArgs;
use crate::{digest_and_filter, kmer_params, marker_source_for_build};
use strain2bscan::cst::SpeciesCst;
use strain2bscan::identify::Params;

pub fn run(args: &DiagnoseTreeArgs) -> Result<(), String> {
    let (k, scale) = kmer_params(Some(args.kmer_size), Some(args.sketch_scale))?;
    let source = marker_source_for_build(&args.marker_source, args.enzyme.as_deref(), k, scale)?;
    let genomes = &args.genomes;
    let similarity = args.similarity;

    let recs = digest_and_filter(genomes, &source, args.max_contigs, args.min_tag_fraction)?;
    let n = recs.len();
    if n < 2 {
        return Err("need at least 2 genomes to form a hierarchy".into());
    }
    let cst = SpeciesCst::build(
        recs.into_iter()
            .map(|r| (r.name, r.markers, r.full_markers))
            .collect(),
        similarity,
        args.containment,
    );
    println!(
        "panel: {} genomes -> {} clusters @ similarity {similarity} ({})",
        n,
        cst.n_clusters(),
        source.describe()
    );

    let stats = cst.hierarchy_stats();
    println!("\n#node\tmembers\tmerge_similarity\tcore\tgroup_specific");
    for s in &stats {
        println!(
            "N{}\t{}\t{:.4}\t{}\t{}",
            s.node_id, s.n_members, s.merge_similarity, s.core, s.group_specific
        );
    }

    let cdb = cst.cluster_db();
    let floor = Params::default().min_support_markers;
    let mut uniq_counts: Vec<usize> = (0..cdb.n_strains())
        .map(|j| cdb.unique_marker_count(j))
        .collect();
    uniq_counts.sort_unstable();
    let below = uniq_counts.iter().filter(|&&u| u < floor).count();
    println!("\n#cluster_unique_markers");
    println!(
        "min={} median={} max={}   below the support floor of {}: {}/{}",
        uniq_counts.first().copied().unwrap_or(0),
        uniq_counts.get(uniq_counts.len() / 2).copied().unwrap_or(0),
        uniq_counts.last().copied().unwrap_or(0),
        floor,
        below,
        uniq_counts.len()
    );
    if below == 0 {
        println!("  -> every cluster already clears the floor on its own markers, so tree pooling");
        println!("     has nothing to rescue here and --layer1 cst will reach the same leaves.");
    } else {
        println!(
            "  -> {below} cluster(s) are invisible to the flat path and only reachable by pooling."
        );
    }

    let mut resolving: Vec<usize> = stats
        .iter()
        .filter(|s| s.merge_similarity < similarity)
        .map(|s| s.group_specific)
        .collect();
    resolving.sort_unstable();
    println!("\n--- verdict ---");
    if resolving.is_empty() {
        println!("no internal node merges below the clustering threshold: the panel is a single cluster,");
        println!("so a tree has nothing to resolve here. Try a more diverse panel.");
        return Ok(());
    }
    let med = resolving[resolving.len() / 2];
    let min = resolving[0];
    println!(
        "{} internal node(s) below the clustering threshold; group-specific markers: min={} median={} max={}",
        resolving.len(),
        min,
        med,
        resolving[resolving.len() - 1]
    );
    const MINK: usize = 25;
    if med >= 4 * MINK {
        println!(
            "VERDICT: tree is viable here (median {med} >= {}). Node sets carry real signal.",
            4 * MINK
        );
    } else if med >= MINK {
        println!(
            "VERDICT: marginal (median {med} in [{MINK}, {})). Viable but with little headroom.",
            4 * MINK
        );
    } else {
        println!("VERDICT: NOT viable (median {med} < {MINK}). Internal nodes are too sparse to");
        println!("descend on; invest in the shared-marker regression instead of the tree.");
    }
    Ok(())
}
