//! Sparse strain × marker database with unique-marker tracking.
//!
//! Unlike `strainscan-rust` (dense `Array2<u8>` serialized to pretty JSON — tens of GB
//! at real scale), this stores, per strain, only the **set of marker hashes it carries**,
//! plus an inverted index `marker -> #strains` so we can flag markers that are unique to
//! a single strain. Unique markers are StrainScan's discriminating signal and, with
//! 2bRAD tags, are exactly the taxonomy-specific tags Fast2bRAD-M already selects in
//! `build_quan_db.rs` (`taxonomies.len() == 1`) — here applied at strain resolution.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::markers::Marker;

#[derive(Debug, Default, Clone)]
pub struct StrainDb {
    pub strain_names: Vec<String>,
    /// Per strain: the set of markers it carries.
    pub strain_markers: Vec<HashSet<Marker>>,
    /// marker -> number of strains carrying it (inverted-index degree).
    pub marker_degree: HashMap<Marker, u32>,
    /// Enzyme set used to build this DB (samples must be digested with the same set).
    pub enzymes: Vec<String>,
    /// Truly-unique markers, defined by **full genome occurrence** (a marker absent — at any
    /// copy number — from every other cluster's genomes), set by CST `cluster_db`. Stricter
    /// than `marker_degree == 1` (single-copy membership), which mislabels a tag as unique when
    /// it is multi-copy in another cluster (single-copy-filter asymmetry) and thus reachable
    /// from that cluster's reads. If empty, uniqueness falls back to `marker_degree`.
    pub unique_set: HashSet<Marker>,
}

impl StrainDb {
    /// Build from `(strain_name, markers)` pairs.
    pub fn build(strains: Vec<(String, Vec<Marker>)>) -> Self {
        let mut db = StrainDb::default();
        for (name, markers) in strains {
            let set: HashSet<Marker> = markers.into_iter().collect();
            for &m in &set {
                *db.marker_degree.entry(m).or_insert(0) += 1;
            }
            db.strain_names.push(name);
            db.strain_markers.push(set);
        }
        db
    }

    pub fn n_strains(&self) -> usize {
        self.strain_names.len()
    }

    /// Is `marker` unique to a single cluster? Uses the occurrence-based `unique_set` when set
    /// (CST databases), else falls back to single-copy membership degree.
    #[inline]
    pub fn is_unique(&self, marker: Marker) -> bool {
        if self.unique_set.is_empty() {
            self.marker_degree.get(&marker).copied() == Some(1)
        } else {
            self.unique_set.contains(&marker)
        }
    }

    /// The unique markers of strain `j`.
    pub fn unique_markers(&self, j: usize) -> impl Iterator<Item = Marker> + '_ {
        self.strain_markers[j]
            .iter()
            .copied()
            .filter(move |&m| self.is_unique(m))
    }

    pub fn unique_marker_count(&self, j: usize) -> usize {
        self.unique_markers(j).count()
    }

    // ===== persistence (simple, line-oriented text) ========================
    // Format:
    //   line 1:            "#strain2bscan-db\t<n_strains>"
    //   next n lines:      "<strain_name>\t<marker_hex,marker_hex,...>"
    // Sparse and compact; production would use a binary/bgzf layout.

    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(
            w,
            "#strain2bscan-db\t{}\t{}",
            self.n_strains(),
            self.enzymes.join(",")
        )?;
        if !self.unique_set.is_empty() {
            let joined = self
                .unique_set
                .iter()
                .map(|m| format!("{m:x}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(w, "#unique\t{joined}")?;
        }
        for (name, markers) in self.strain_names.iter().zip(&self.strain_markers) {
            let joined = markers
                .iter()
                .map(|m| format!("{m:x}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(w, "{name}\t{joined}")?;
        }
        Ok(())
    }

    pub fn load(path: &Path) -> std::io::Result<Self> {
        let reader = BufReader::new(File::open(path)?);
        let mut strains = Vec::new();
        let mut enzymes: Vec<String> = Vec::new();
        let mut unique_set: HashSet<Marker> = HashSet::new();
        let hexset = |csv: &str| {
            csv.split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| Marker::from_str_radix(s, 16).ok())
                .collect::<HashSet<Marker>>()
        };
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                if line.starts_with("#strain2bscan-db") {
                    if let Some(csv) = line.split('\t').nth(2) {
                        enzymes = csv.split(',').filter(|s| !s.is_empty()).map(String::from).collect();
                    }
                } else if line.starts_with("#unique") {
                    if let Some(csv) = line.split('\t').nth(1) {
                        unique_set = hexset(csv);
                    }
                }
                continue;
            }
            if line.is_empty() {
                continue;
            }
            let mut it = line.splitn(2, '\t');
            let name = it.next().unwrap_or("").to_string();
            let markers = it
                .next()
                .unwrap_or("")
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| Marker::from_str_radix(s, 16).ok())
                .collect::<Vec<_>>();
            strains.push((name, markers));
        }
        let mut db = StrainDb::build(strains);
        db.enzymes = enzymes;
        db.unique_set = unique_set;
        Ok(db)
    }

    /// Quick DB stats for the `info`/`build` CLI.
    pub fn stats(&self) -> DbStats {
        let total_markers = self.marker_degree.len();
        let unique_total = self.marker_degree.values().filter(|&&d| d == 1).count();
        let avg_markers = if self.n_strains() == 0 {
            0.0
        } else {
            self.strain_markers.iter().map(|s| s.len()).sum::<usize>() as f64
                / self.n_strains() as f64
        };
        DbStats {
            n_strains: self.n_strains(),
            n_markers: total_markers,
            unique_markers: unique_total,
            avg_markers_per_strain: avg_markers,
            unique_fraction: if total_markers == 0 {
                0.0
            } else {
                unique_total as f64 / total_markers as f64
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct DbStats {
    pub n_strains: usize,
    pub n_markers: usize,
    pub unique_markers: usize,
    pub avg_markers_per_strain: f64,
    pub unique_fraction: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> StrainDb {
        // markers 1,2,3 shared "core"; 10/20/30 are private to A/B/C.
        StrainDb::build(vec![
            ("A".into(), vec![1, 2, 3, 10]),
            ("B".into(), vec![1, 2, 3, 20]),
            ("C".into(), vec![1, 2, 3, 30]),
        ])
    }

    #[test]
    fn unique_markers_are_identified() {
        let db = toy();
        assert!(db.is_unique(10) && db.is_unique(20) && db.is_unique(30));
        assert!(!db.is_unique(1));
        assert_eq!(db.unique_marker_count(0), 1);
        assert_eq!(db.unique_markers(0).next(), Some(10));
    }

    #[test]
    fn roundtrip_save_load() {
        let db = toy();
        let dir = std::env::temp_dir();
        let path = dir.join("strain2bscan_test_db.tsv");
        db.save(&path).unwrap();
        let back = StrainDb::load(&path).unwrap();
        assert_eq!(back.n_strains(), 3);
        assert!(back.is_unique(20));
        let _ = std::fs::remove_file(path);
    }
}
