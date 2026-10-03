//! Multi-species panel loading, species gating, and per-sample profiling.

use std::path::{Path, PathBuf};

use strain2bscan::db::StrainDb;
use strain2bscan::enzymes::parse_enzyme_set;
use strain2bscan::fxhash::{FxHashMap, FxHashSet};
use strain2bscan::identify::{
    detectable_fraction, min_count_for, profile, resolve_layer1, tree_utility, Layer1, Params,
    StrainCall, TreeUtility,
};
use strain2bscan::markers::{parse_kmer_db_token, Marker, MarkerCounts};
use strain2bscan::parallel::par_map;

use crate::MarkerSource;

/// The three-tier species-gate thresholds.
#[derive(Clone, Copy)]
pub struct SpeciesGate {
    pub min_markers: usize,
    pub min_frac: f64,
    pub min_detect: usize,
}

/// Layer-1 outcome for one species in a multi-species sample.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SpeciesTier {
    /// Enough species-specific marker evidence to attempt strain resolution.
    Resolved,
    /// Present, but below the evidence needed to resolve strains.
    DetectedNotResolved,
    /// Not enough evidence to call the species present.
    Absent,
}

/// Classify a species from ABSOLUTE species-specific marker evidence.
///
/// `present` = species-specific markers observed; `total` = species-specific markers in the DB;
/// `detect`/`floor` are absolute thresholds and `frac` is the breadth fraction; `reachable` is
/// the fraction of the panel observable at this species' estimated depth.
pub fn species_tier(
    present: usize,
    total: usize,
    detect: usize,
    floor: usize,
    frac: f64,
    reachable: f64,
) -> SpeciesTier {
    const MIN_FLOOR_FRACTION: f64 = 0.25;
    let r = reachable.clamp(0.0, 1.0).max(MIN_FLOOR_FRACTION);
    let frac_gate = (frac.max(0.0) * total as f64 * r).ceil() as usize;
    let scaled_floor = (floor as f64 * r).ceil() as usize;
    let resolve_gate = scaled_floor.max(frac_gate).max(detect).max(1);
    let detect_gate = detect.min(resolve_gate);
    if present >= resolve_gate {
        SpeciesTier::Resolved
    } else if present >= detect_gate {
        SpeciesTier::DetectedNotResolved
    } else {
        SpeciesTier::Absent
    }
}

pub fn species_gate(
    min_markers: Option<usize>,
    min_frac: Option<f64>,
    min_detect: Option<usize>,
) -> SpeciesGate {
    SpeciesGate {
        min_markers: min_markers.unwrap_or(200),
        min_frac: min_frac.unwrap_or(0.0),
        min_detect: min_detect.unwrap_or(10),
    }
}

/// Per-species Layer-1 result carried out of the parallel map.
pub struct SpeciesResult {
    pub species: String,
    pub present_specific: usize,
    pub total_specific: usize,
    /// Estimated per-tag depth over this species' species-specific markers.
    pub lambda: f64,
    pub tier: SpeciesTier,
    pub layer1: Layer1,
    pub tree: Option<TreeUtility>,
    pub calls: Vec<StrainCall>,
}

/// A loaded multi-species panel.
pub struct Panel {
    pub loaded: Vec<(String, StrainDb)>,
    pub source: MarkerSource,
    pub specific_sets: Vec<FxHashSet<Marker>>,
}

/// Collect + load the per-species DBs, decide the sample marker source, and apply the
/// cross-species marker restriction.
pub fn load_panel(
    dbs_dir: &Path,
    marker_source_arg: Option<&str>,
    enzyme_arg: Option<&str>,
    kmer_size_arg: Option<usize>,
    sketch_scale_arg: Option<u64>,
    no_cross_species_filter: bool,
) -> Result<Panel, String> {
    let mut loaded = load_species_dbs(dbs_dir)?;

    let kmer_of = |db: &StrainDb| {
        if db.enzymes.len() == 1 {
            parse_kmer_db_token(&db.enzymes[0])
        } else {
            None
        }
    };
    let n_kmer = loaded
        .iter()
        .filter(|(_, db)| kmer_of(db).is_some())
        .count();

    let source: MarkerSource = if n_kmer == 0 {
        if marker_source_arg == Some("kmer") {
            return Err(
                "--marker-source kmer but every DB in --dbs is an enzyme-tag database; \
                 the marker spaces are disjoint"
                    .into(),
            );
        }
        let Some(spec) = enzyme_arg else {
            return Err("missing --enzyme".into());
        };
        let set = parse_enzyme_set(spec).ok_or_else(|| format!("unknown enzyme set: {spec}"))?;
        MarkerSource::Enzyme(set)
    } else if n_kmer != loaded.len() {
        return Err(format!(
            "mixed panel: {} k-mer DB(s) and {} enzyme-tag DB(s) in --dbs; a sample can only \
             be digested one way per run",
            n_kmer,
            loaded.len() - n_kmer
        ));
    } else {
        let (k, scale) = kmer_of(&loaded[0].1).unwrap();
        if let Some((sp, _)) = loaded
            .iter()
            .find(|(_, db)| kmer_of(db) != Some((k, scale)))
        {
            return Err(format!(
                "k-mer DBs in --dbs disagree on k/scale: '{sp}' was not built with \
                 k={k}, scale={scale}"
            ));
        }
        if marker_source_arg == Some("enzyme") {
            return Err(format!(
                "this panel is a k-mer sketch (k={k}, scale={scale}); --marker-source enzyme \
                 would compare disjoint marker spaces"
            ));
        }
        if let Some(want_k) = kmer_size_arg {
            if want_k != k {
                return Err(format!(
                    "--kmer-size {want_k} does not match the panel (built with k={k})"
                ));
            }
        }
        if let Some(want_scale) = sketch_scale_arg {
            if want_scale != scale {
                return Err(format!(
                    "--sketch-scale {want_scale} does not match the panel (built with scale={scale})"
                ));
            }
        }
        if let Some(e) = enzyme_arg {
            eprintln!("warning: --enzyme {e} is ignored: the panel is a k-mer sketch (k={k}, scale={scale})");
        }
        MarkerSource::Kmer { k, scale }
    };

    let mut species_degree: FxHashMap<Marker, u32> = FxHashMap::default();
    for (_, db) in &loaded {
        for &m in db.marker_degree.keys() {
            *species_degree.entry(m).or_insert(0) += 1;
        }
    }

    let specific_sets: Vec<FxHashSet<Marker>> = loaded
        .iter()
        .map(|(_, db)| {
            db.marker_degree
                .keys()
                .copied()
                .filter(|m| species_degree.get(m).copied() == Some(1))
                .collect()
        })
        .collect();

    if !no_cross_species_filter {
        let (mut before, mut after) = (0usize, 0usize);
        for ((_, db), specific) in loaded.iter_mut().zip(&specific_sets) {
            before += db.marker_degree.len();
            db.restrict_to(specific);
            after += specific.len();
        }
        println!(
            "cross-species filter: {after}/{before} markers usable for quantification ({:.1}% shared with another species in the panel, excluded)",
            if before > 0 {
                100.0 * (before - after) as f64 / before as f64
            } else {
                0.0
            }
        );
    }

    Ok(Panel {
        loaded,
        source,
        specific_sets,
    })
}

/// Collect and load every per-species DB in `dbs_dir`.
pub fn load_species_dbs(dbs_dir: &Path) -> Result<Vec<(String, StrainDb)>, String> {
    let mut db_paths: Vec<PathBuf> = std::fs::read_dir(dbs_dir)
        .map_err(|e| format!("cannot list DB dir {}: {e}", dbs_dir.display()))?
        .map(|e| {
            e.map(|e| e.path()).map_err(|err| {
                format!(
                    "cannot read an entry of DB dir {}: {err}",
                    dbs_dir.display()
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|p| {
            p.extension()
                .and_then(|x| x.to_str())
                .is_some_and(|x| x.eq_ignore_ascii_case("tsv"))
                && p.file_name()
                    .and_then(|x| x.to_str())
                    .is_some_and(|n| !n.contains(".members."))
        })
        .collect();
    db_paths.sort();
    if db_paths.is_empty() {
        return Err("no *.tsv species DBs found in --dbs dir".into());
    }
    par_map(&db_paths, |path| {
        let sp = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        StrainDb::load(path)
            .map(|db| (sp, db))
            .map_err(|e| format!("failed to load species DB {}: {e}", path.display()))
    })
    .into_iter()
    .collect()
}

/// One output row of a multi-species profile.
pub struct ProfileRow {
    pub species: String,
    pub call: StrainCall,
    pub global_abundance: f64,
    pub sample_fraction: f64,
}

/// Everything one sample yields against the panel, before any reporting.
pub struct SampleProfile {
    pub per_species: Vec<SpeciesResult>,
    pub rows: Vec<ProfileRow>,
    pub total_tags: u64,
}

/// Gate + strain-profile every species of the panel against ONE sample's tag counts.
pub fn profile_sample(
    panel: &Panel,
    counts: &MarkerCounts,
    gate: &SpeciesGate,
    params: &Params,
    min_global_abundance: Option<f64>,
) -> Result<SampleProfile, String> {
    let loaded = &panel.loaded;
    let specific_sets = &panel.specific_sets;
    let gate = *gate;

    let order: Vec<usize> = (0..loaded.len()).collect();
    let mut per_species: Vec<SpeciesResult> = par_map(&order, |&i| {
        let (species, db) = &loaded[i];
        let specific = &specific_sets[i];
        let total_specific = specific.len();
        let observed: u64 = specific
            .iter()
            .map(|m| counts.get(m).copied().unwrap_or(0) as u64)
            .sum();
        let lambda = if total_specific > 0 {
            observed as f64 / total_specific as f64
        } else {
            0.0
        };
        let min_count = if params.adaptive_singleton {
            min_count_for(lambda)
        } else {
            2
        };
        let present_specific = specific
            .iter()
            .filter(|m| counts.get(*m).copied().unwrap_or(0) >= min_count)
            .count();
        let reachable = if params.adaptive_floor {
            detectable_fraction(lambda)
        } else {
            1.0
        };
        let tier = species_tier(
            present_specific,
            total_specific,
            gate.min_detect,
            gate.min_markers,
            gate.min_frac,
            reachable,
        );
        let calls = if tier == SpeciesTier::Resolved {
            profile(db, counts, params)
        } else {
            Vec::new()
        };
        let layer1 = resolve_layer1(db, params);
        let tree = tree_utility(db, params.min_support_markers);
        SpeciesResult {
            species: species.clone(),
            present_specific,
            total_specific,
            lambda,
            layer1,
            tree,
            tier,
            calls,
        }
    });

    if let Some(min_global) = min_global_abundance {
        if min_global > 0.0 {
            let total: f64 = per_species
                .iter()
                .flat_map(|r| r.calls.iter())
                .map(|c| c.depth)
                .sum();
            if total > 0.0 {
                for r in &mut per_species {
                    r.calls.retain(|c| c.depth / total >= min_global);
                }
                for r in &mut per_species {
                    if r.tier == SpeciesTier::Resolved && r.calls.is_empty() {
                        r.tier = SpeciesTier::DetectedNotResolved;
                    }
                }
                for r in &mut per_species {
                    let kept: f64 = r.calls.iter().map(|c| c.rel_abundance).sum();
                    if kept > 0.0 {
                        for c in &mut r.calls {
                            c.rel_abundance /= kept;
                        }
                    }
                }
            }
        }
    }

    let depth_sum: f64 = per_species
        .iter()
        .flat_map(|r| r.calls.iter())
        .map(|c| c.depth)
        .sum();
    let global_of = |c: &StrainCall| {
        if depth_sum > 0.0 {
            c.depth / depth_sum
        } else {
            0.0
        }
    };

    let total_tags: u64 = counts.values().map(|&c| c as u64).sum();
    let mass_of = |c: &StrainCall| c.depth * c.n_markers as f64;
    let sample_fraction_of = |c: &StrainCall| {
        if total_tags > 0 {
            mass_of(c) / total_tags as f64
        } else {
            0.0
        }
    };

    let mut rows: Vec<ProfileRow> = per_species
        .iter()
        .flat_map(|r| {
            r.calls.iter().map(move |c| ProfileRow {
                species: r.species.clone(),
                call: c.clone(),
                global_abundance: global_of(c),
                sample_fraction: sample_fraction_of(c),
            })
        })
        .collect();
    rows.sort_by(|a, b| {
        a.species.cmp(&b.species).then(
            b.call
                .rel_abundance
                .partial_cmp(&a.call.rel_abundance)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
    });

    Ok(SampleProfile {
        per_species,
        rows,
        total_tags,
    })
}

#[cfg(test)]
mod tests {
    use super::{load_species_dbs, species_tier, SpeciesTier};
    use strain2bscan::db::StrainDb;
    use strain2bscan::depth::detectable_fraction;

    fn kmer_db() -> StrainDb {
        let mut db = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
        ]);
        db.enzymes = vec!["kmer15s1".to_string()];
        db
    }

    fn enzyme_db() -> StrainDb {
        let mut db = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
        ]);
        db.enzymes = vec!["BcgI".to_string()];
        db
    }

    #[test]
    fn multi_profile_panel_load_failure_propagates() {
        let dir = std::env::temp_dir().join(format!("s2bs_panel_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let good = StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
        ]);
        good.save(&dir.join("good.tsv")).unwrap();
        std::fs::write(
            dir.join("broken.tsv"),
            b"#strain2bscan-db\t5\t\t1,1,1,1,1\nX\t1,2,3\n",
        )
        .unwrap();

        let err = load_species_dbs(&dir).unwrap_err();
        assert!(
            err.contains("broken.tsv"),
            "error must name the failing DB: {err}"
        );

        std::fs::remove_file(dir.join("broken.tsv")).unwrap();
        let loaded = load_species_dbs(&dir).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].0, "good");
        assert_eq!(loaded[0].1.n_strains(), 2);
        let _ = std::fs::remove_dir_all(&dir);
    }

    const FULL: f64 = 1.0;

    #[test]
    fn absolute_floor_gates_when_no_fraction() {
        assert_eq!(
            species_tier(250, 5000, 10, 200, 0.0, FULL),
            SpeciesTier::Resolved
        );
        assert_eq!(
            species_tier(50, 5000, 10, 200, 0.0, FULL),
            SpeciesTier::DetectedNotResolved
        );
        assert_eq!(
            species_tier(5, 5000, 10, 200, 0.0, FULL),
            SpeciesTier::Absent
        );
    }

    #[test]
    fn breadth_fraction_raises_the_bar_for_large_panels() {
        assert_eq!(
            species_tier(300, 5000, 10, 200, 0.10, FULL),
            SpeciesTier::DetectedNotResolved
        );
        assert_eq!(
            species_tier(600, 5000, 10, 200, 0.10, FULL),
            SpeciesTier::Resolved
        );
    }

    #[test]
    fn small_panel_species_can_still_be_detected() {
        assert_eq!(
            species_tier(150, 150, 10, 200, 0.0, FULL),
            SpeciesTier::DetectedNotResolved
        );
        assert_eq!(
            species_tier(5, 150, 10, 200, 0.0, FULL),
            SpeciesTier::Absent
        );
    }

    #[test]
    fn low_depth_relaxes_the_floor_but_only_to_the_bound() {
        let reachable = detectable_fraction(0.05);
        assert!(reachable < 0.05);
        assert_eq!(
            species_tier(60, 5000, 10, 200, 0.0, FULL),
            SpeciesTier::DetectedNotResolved
        );
        assert_eq!(
            species_tier(60, 5000, 10, 200, 0.0, reachable),
            SpeciesTier::Resolved
        );
    }

    #[test]
    fn adaptive_floor_does_not_cancel_itself_away() {
        let r = detectable_fraction(0.002);
        assert_eq!(
            species_tier(40, 20_000, 10, 200, 0.0, r),
            SpeciesTier::DetectedNotResolved
        );
        assert_eq!(
            species_tier(50, 20_000, 10, 200, 0.0, r),
            SpeciesTier::Resolved
        );
    }

    #[test]
    fn adaptive_gate_never_falls_below_the_detect_floor() {
        for lambda in [0.0, 1e-6, 0.001, 0.01] {
            let r = detectable_fraction(lambda);
            assert_eq!(species_tier(9, 5000, 10, 200, 0.0, r), SpeciesTier::Absent);
            assert_eq!(
                species_tier(10, 5000, 10, 200, 0.0, r),
                SpeciesTier::DetectedNotResolved
            );
        }
    }
}
