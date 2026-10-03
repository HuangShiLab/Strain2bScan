//! Cluster Search Tree demo.

use std::collections::HashMap;

use strain2bscan::cst::{SpeciesCst, DEFAULT_SIMILARITY};
use strain2bscan::identify::{profile, Params};
use strain2bscan::markers::{Marker, MarkerCounts};

use crate::report::report;

pub fn run() -> Result<(), String> {
    let core: Vec<Marker> = (0..200).collect();
    let clu_a: Vec<Marker> = (200..240).collect();
    let clu_b: Vec<Marker> = (300..340).collect();
    let mk = |name: &str, extra: &[Marker], base: Marker| {
        let mut v = core.clone();
        v.extend_from_slice(extra);
        v.extend((0..3).map(|i| base + i));
        (name.to_string(), v.clone(), v)
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
