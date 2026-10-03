//! Build a Cluster Search Tree DB from genomes.

use crate::cli::ClusterArgs;
use crate::{digest_and_filter, kmer_params, marker_source_for_build};
use strain2bscan::cst::{SpeciesCst, MINHASH_ABOVE};
use strain2bscan::identify::Params;
use strain2bscan::parallel::num_threads;

pub fn run(args: &ClusterArgs) -> Result<(), String> {
    let (k, scale) = kmer_params(Some(args.kmer_size), Some(args.sketch_scale))?;
    let source = marker_source_for_build(&args.marker_source, args.enzyme.as_deref(), k, scale)?;
    let genomes = &args.genomes;
    let out = &args.out;
    let similarity = args.similarity;
    let containment = args.containment;

    let recs = digest_and_filter(genomes, &source, args.max_contigs, args.min_tag_fraction)?;
    let n_genomes = recs.len();
    let cst = SpeciesCst::build(
        recs.into_iter()
            .map(|r| (r.name, r.markers, r.full_markers))
            .collect(),
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
        "clustered {} genomes into {} cluster(s) @ similarity {similarity} ({}, threads: {}, clustering: {}-{})",
        cst.genome_names.len(),
        cst.n_clusters(),
        source.describe(),
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

    let mut db = cst.cluster_db();
    db.enzymes = source.db_token();
    db.tree = Some(cst.build_tree());
    let min_markers = Params::default().min_support_markers;
    let mut resolvable = 0usize;
    let with_what = match &source {
        crate::MarkerSource::Enzyme(_) => "with this enzyme set".to_string(),
        crate::MarkerSource::Kmer { k, scale } => {
            format!("with this k-mer sketch (k={k}, scale={scale})")
        }
    };
    for cid in 0..db.n_strains() {
        let n_spec = db.unique_marker_count(cid);
        if n_spec >= min_markers {
            resolvable += 1;
        } else {
            println!(
                "  ⚠ C{cid} has only {n_spec} cluster-specific markers (< {min_markers}); \
                 not reliably resolvable {with_what}."
            );
        }
    }
    if resolvable == 0 {
        match &source {
            crate::MarkerSource::Enzyme(set) => println!(
                "  ✗ NOT DOABLE at strain/cluster level for this species with enzyme(s) {}. \
                 The species can still be detected (Layer-1); for finer resolution use more \
                 enzymes (--enzyme all) on a conventional metagenome.",
                set.iter().map(|e| e.name).collect::<Vec<_>>().join("+")
            ),
            crate::MarkerSource::Kmer { .. } => println!(
                "  ✗ NOT DOABLE at strain/cluster level for this species {with_what}. \
                 The species can still be detected (Layer-1); for finer resolution use a \
                 smaller k or a denser sketch (lower --sketch-scale)."
            ),
        }
    }

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

    db.save(out).map_err(|e| e.to_string())?;
    println!(
        "saved cluster DB -> {} ({} clusters, {} resolvable); membership -> {}",
        out.display(),
        db.n_strains(),
        resolvable,
        members_path.display()
    );
    Ok(())
}
