//! In-memory conspecific demo.

use std::collections::HashMap;

use strain2bscan::db::StrainDb;
use strain2bscan::identify::{naive_profile, profile, Params};
use strain2bscan::markers::{Marker, MarkerCounts};

use crate::report::report;

pub fn run() -> Result<(), String> {
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
