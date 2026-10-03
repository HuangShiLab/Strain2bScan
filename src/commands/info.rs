//! Print information about a strain DB.

use crate::cli::InfoArgs;
use crate::print_stats;
use strain2bscan::db::StrainDb;
use strain2bscan::markers::parse_kmer_db_token;

pub fn run(args: &InfoArgs) -> Result<(), String> {
    let db = StrainDb::load(&args.db).map_err(|e| e.to_string())?;
    let kmer = if db.enzymes.len() == 1 {
        parse_kmer_db_token(&db.enzymes[0])
    } else {
        None
    };
    match kmer {
        Some((k, scale)) => println!(
            "marker source: k-mer sketch (k={k}, sketch scale={scale}; DB token: {})",
            db.enzymes[0]
        ),
        None => println!(
            "enzymes: {}",
            if db.enzymes.is_empty() {
                "(unspecified)".into()
            } else {
                db.enzymes.join("+")
            }
        ),
    }
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
