//! Markers = canonical 2bRAD tags hashed to `u64`.
//!
//! A marker is the hash of the canonical (lexicographically smaller of forward /
//! reverse-complement) tag sequence. Hashing keeps the representation length-agnostic
//! (tags are 32–38 bp depending on enzyme), which is why we use `u64` markers rather
//! than `strainscan-rust`'s 2-bit-packed 31-mers (those cap at 31 bp / need u128 here).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::enzymes::Enzyme;

/// A strain/sample marker: hash of a canonical 2bRAD tag.
pub type Marker = u64;

/// Reverse-complement, preserving N.
#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

/// Canonical orientation: min(seq, revcomp(seq)) as raw bytes (matches Fast2bRAD-M).
pub fn canonical(seq: &[u8]) -> Vec<u8> {
    let rc = revcomp(seq);
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

/// FNV-1a 64-bit hash.
///
/// NOTE: Fast2bRAD-M hashes tags with `fxhash::FxHasher`. To consume its `*.iibdb`
/// files directly, production should hash with the identical FxHash function instead
/// of FNV so marker values are byte-compatible. We use FNV here only to stay
/// dependency-free; DBs built by this prototype are internally consistent.
pub fn hash_bytes(bytes: &[u8]) -> u64 {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut h = OFFSET;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(PRIME);
    }
    h
}

/// Tag sequence (any orientation/case) → canonical marker.
pub fn marker_from_tag(tag: &[u8]) -> Marker {
    let mut up = tag.to_vec();
    up.make_ascii_uppercase();
    hash_bytes(&canonical(&up))
}

/// Digest one sequence, returning every tag marker (with multiplicity).
pub fn digest_sequence(seq: &[u8], enzyme: &Enzyme) -> Vec<Marker> {
    let mut up = seq.to_vec();
    up.make_ascii_uppercase();
    enzyme
        .find_all_tags(&up)
        .into_iter()
        .map(|(pos, len)| marker_from_tag(&up[pos..pos + len]))
        .collect()
}

/// Digest one sequence with a **set** of enzymes, pooling all tag markers.
///
/// Used for conventional metagenomes: digitally digesting reads with all 16 type-IIB
/// enzymes enriches the marker pool ~16× vs. BcgI alone, recovering more strain-specific
/// loci. Different enzymes yield different-length tags → distinct hashes, so pooling is safe.
pub fn digest_sequence_multi(seq: &[u8], enzymes: &[&Enzyme]) -> Vec<Marker> {
    let mut up = seq.to_vec();
    up.make_ascii_uppercase();
    let mut out = Vec::new();
    for enzyme in enzymes {
        for (pos, len) in enzyme.find_all_tags(&up) {
            out.push(marker_from_tag(&up[pos..pos + len]));
        }
    }
    out
}

// ===== minimal FASTA / FASTQ reading (prototype only) ======================
// Production: replace with needletail (streaming + gzip). These load whole files.

/// Read sequences from a FASTA file (concatenating wrapped lines per record).
pub fn read_fasta(path: &Path) -> std::io::Result<Vec<Vec<u8>>> {
    let reader = BufReader::new(File::open(path)?);
    let mut seqs = Vec::new();
    let mut cur = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !cur.is_empty() {
                seqs.push(std::mem::take(&mut cur));
            }
        } else {
            cur.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    Ok(seqs)
}

/// Read sequences from a FASTQ file (line 2 of every 4-line record).
pub fn read_fastq(path: &Path) -> std::io::Result<Vec<Vec<u8>>> {
    let reader = BufReader::new(File::open(path)?);
    let mut seqs = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        if i % 4 == 1 {
            seqs.push(line?.trim().as_bytes().to_vec());
        }
    }
    Ok(seqs)
}

/// Read FASTA or FASTQ by extension.
pub fn read_fastx(path: &Path) -> std::io::Result<Vec<Vec<u8>>> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("fq") | Some("fastq") => read_fastq(path),
        _ => read_fasta(path),
    }
}

/// Digest all sequences and return the **set** of markers (for a reference genome).
pub fn genome_markers(seqs: &[Vec<u8>], enzyme: &Enzyme) -> Vec<Marker> {
    let mut set = std::collections::HashSet::new();
    for s in seqs {
        for m in digest_sequence(s, enzyme) {
            set.insert(m);
        }
    }
    set.into_iter().collect()
}

/// Per-genome marker copy numbers (how many times each tag occurs in the genome).
pub fn genome_marker_counts(seqs: &[Vec<u8>], enzyme: &Enzyme) -> HashMap<Marker, u32> {
    let mut counts: HashMap<Marker, u32> = HashMap::new();
    for s in seqs {
        for m in digest_sequence(s, enzyme) {
            *counts.entry(m).or_insert(0) += 1;
        }
    }
    counts
}

/// Single-copy markers only (tags occurring exactly once in the genome).
///
/// StrainScan and Fast2bRAD-M (`remove_redundant`) both restrict markers to single-copy
/// loci, because multi-copy tags inflate abundance estimates and blur strain identity.
pub fn single_copy_markers(counts: &HashMap<Marker, u32>) -> Vec<Marker> {
    counts
        .iter()
        .filter(|(_, &c)| c == 1)
        .map(|(&m, _)| m)
        .collect()
}

/// Digest all reads and return marker **counts** (for a sample).
pub fn sample_marker_counts(seqs: &[Vec<u8>], enzyme: &Enzyme) -> HashMap<Marker, u32> {
    sample_marker_counts_multi(seqs, &[enzyme])
}

/// Per-genome marker copy numbers with a set of enzymes (pooled).
pub fn genome_marker_counts_multi(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> HashMap<Marker, u32> {
    let mut counts: HashMap<Marker, u32> = HashMap::new();
    for s in seqs {
        for m in digest_sequence_multi(s, enzymes) {
            *counts.entry(m).or_insert(0) += 1;
        }
    }
    counts
}

/// Sample marker counts with a set of enzymes (pooled).
pub fn sample_marker_counts_multi(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> HashMap<Marker, u32> {
    genome_marker_counts_multi(seqs, enzymes)
}

/// Parallel sample marker counts: digest read chunks across threads, then merge the maps.
/// Read digestion dominates per-sample profiling time, so this is the main speedup path.
pub fn sample_marker_counts_multi_par(
    seqs: &[Vec<u8>],
    enzymes: &[&Enzyme],
) -> HashMap<Marker, u32> {
    let nt = crate::parallel::num_threads().min(seqs.len().max(1));
    if nt <= 1 || seqs.len() < 4096 {
        return genome_marker_counts_multi(seqs, enzymes);
    }
    let chunk = seqs.len().div_ceil(nt);
    let chunks: Vec<&[Vec<u8>]> = seqs.chunks(chunk).collect();
    let partials: Vec<HashMap<Marker, u32>> =
        crate::parallel::par_map(&chunks, |c| genome_marker_counts_multi(c, enzymes));
    let mut out: HashMap<Marker, u32> = HashMap::new();
    for p in partials {
        for (m, c) in p {
            *out.entry(m).or_insert(0) += c;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_is_strand_invariant() {
        let fwd = b"ACGTTGCA";
        let rev = revcomp(fwd);
        assert_eq!(canonical(fwd), canonical(&rev));
        assert_eq!(marker_from_tag(fwd), marker_from_tag(&rev));
    }

    #[test]
    fn hash_is_deterministic_and_case_insensitive() {
        assert_eq!(marker_from_tag(b"acgtACGT"), marker_from_tag(b"ACGTACGT"));
        assert_ne!(marker_from_tag(b"AAAA"), marker_from_tag(b"AAAC"));
    }

    #[test]
    fn digest_counts_multiplicity() {
        use crate::enzymes::BCGI;
        // Two identical BcgI sites → the same canonical marker counted twice.
        let mut w = vec![b'A'; 32];
        w[10..13].copy_from_slice(b"CGA");
        w[19..22].copy_from_slice(b"TGC");
        let mut seq = w.clone();
        seq.extend_from_slice(b"GG");
        seq.extend_from_slice(&w);
        let counts = sample_marker_counts(&[seq], &BCGI);
        assert!(counts.values().any(|&c| c >= 2));
    }
}
