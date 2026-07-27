//! Markers = canonical 2bRAD tags hashed to `u64`.
//!
//! A marker is the hash of the canonical (lexicographically smaller of forward /
//! reverse-complement) tag sequence. Hashing keeps the representation length-agnostic
//! (tags are 25–33 bp depending on enzyme), which is why we use `u64` markers rather
//! than `strainscan-rust`'s 2-bit-packed 31-mers (those cap at 31 bp / need u128 here).
//!
//! ## Hot path
//!
//! Digestion is the dominant cost of both DB build and profiling, so the path
//! sequence → tag → marker → count is allocation-free:
//!   * enzyme scanning is case-insensitive ([`Enzyme::for_each_tag`]), so sequences are read
//!     straight from the input buffer — no upper-cased copy per read/contig;
//!   * canonicalization ([`marker_from_tag`]) decides orientation by streaming comparison and
//!     hashes the winner in place — no `revcomp` buffer, no `to_vec`;
//!   * counts land directly in an [`FxHashMap`] — no intermediate `Vec<Marker>` per sequence.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::process::{Child, Command, Stdio};

use crate::enzymes::Enzyme;
use crate::fxhash::{FxHashMap, FxHashSet};

/// A strain/sample marker: hash of a canonical 2bRAD tag.
pub type Marker = u64;

/// Marker → observed count. Keyed by `u64`, so it uses FxHash rather than SipHash.
pub type MarkerCounts = FxHashMap<Marker, u32>;

/// Longest tag any enzyme in the table emits is 33 bp (CspCI); 48 leaves headroom and keeps
/// the stack buffer within one cache line pair.
pub const MAX_TAG_LEN: usize = 48;

/// Complement one base, preserving anything that is not ACGT (e.g. N).
#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

/// Reverse-complement, preserving N.
#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp(b)).collect()
}

/// Canonical orientation: min(seq, revcomp(seq)) as raw bytes (matches Fast2bRAD-M).
///
/// Kept as the readable reference definition; the hot path uses [`marker_from_tag`], which
/// computes the same value without materializing either strand.
pub fn canonical(seq: &[u8]) -> Vec<u8> {
    let rc = revcomp(seq);
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

/// FNV-1a 64-bit hash. Dependency-free; genome and sample tags are hashed with the same
/// function, so marker values are internally consistent across DB build and profiling.
///
/// This defines the on-disk marker identity — do not change it without invalidating every
/// database built by an earlier version.
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

/// FNV-1a of `revcomp(seq)`, computed without building the reverse complement.
#[inline]
fn hash_bytes_rc(seq: &[u8]) -> u64 {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut h = OFFSET;
    for &b in seq.iter().rev() {
        h ^= comp(b) as u64;
        h = h.wrapping_mul(PRIME);
    }
    h
}

/// Tag sequence (any orientation/case) → canonical marker. Allocation-free.
///
/// Equivalent to `hash_bytes(&canonical(&uppercase(tag)))`, but the orientation is chosen by
/// comparing the forward strand against its reverse complement one base at a time (the RC of
/// position `i` is `comp(tag[n-1-i])`), then the winning strand is hashed in place.
#[inline]
pub fn marker_from_tag(tag: &[u8]) -> Marker {
    let n = tag.len();
    if n > MAX_TAG_LEN {
        // Not reachable from the enzyme table; kept correct for external callers.
        let mut up = tag.to_vec();
        up.make_ascii_uppercase();
        return hash_bytes(&canonical(&up));
    }
    let mut buf = [0u8; MAX_TAG_LEN];
    for (i, &b) in tag.iter().enumerate() {
        buf[i] = b.to_ascii_uppercase();
    }
    let up = &buf[..n];

    // Decide min(up, revcomp(up)) without materializing the reverse complement.
    let mut forward_wins = true;
    for i in 0..n {
        let a = up[i];
        let b = comp(up[n - 1 - i]);
        if a != b {
            forward_wins = a < b;
            break;
        }
    }
    if forward_wins {
        hash_bytes(up)
    } else {
        hash_bytes_rc(up)
    }
}

/// Digest one sequence, returning every tag marker (with multiplicity).
///
/// Forward-strand scan only. Each enzyme carries forward+reverse patterns that are exact
/// reverse-complement pairs at its tag length, so a single forward pass finds sites in both
/// orientations, and `marker_from_tag` canonicalizes each tag — making digestion strand-invariant
/// (`digest(seq) == digest(revcomp(seq))`) without scanning both strands. This mirrors
/// Fast2bRAD-M / `2bRADExtraction.pl` and keeps the marker set at 1× (no doubling).
pub fn digest_sequence(seq: &[u8], enzyme: &Enzyme) -> Vec<Marker> {
    let mut out = Vec::new();
    enzyme.for_each_tag(seq, |pos, len| out.push(marker_from_tag(&seq[pos..pos + len])));
    out
}

/// Digest one sequence with a **set** of enzymes, pooling all tag markers.
///
/// Used for conventional metagenomes: digitally digesting reads with all 16 type-IIB
/// enzymes enriches the marker pool ~16× vs. BcgI alone, recovering more strain-specific
/// loci. Different enzymes yield different-length tags → distinct hashes, so pooling is safe.
pub fn digest_sequence_multi(seq: &[u8], enzymes: &[&Enzyme]) -> Vec<Marker> {
    let mut out = Vec::new();
    for enzyme in enzymes {
        enzyme.for_each_tag(seq, |pos, len| out.push(marker_from_tag(&seq[pos..pos + len])));
    }
    out
}

/// Digest one sequence into an existing count map. The allocation-free hot path.
#[inline]
pub fn count_markers_into(seq: &[u8], enzymes: &[&Enzyme], counts: &mut MarkerCounts) {
    for enzyme in enzymes {
        enzyme.for_each_tag(seq, |pos, len| {
            *counts.entry(marker_from_tag(&seq[pos..pos + len])).or_insert(0) += 1;
        });
    }
}

// ===== streaming FASTA / FASTQ reading (plain or gzip) =====================

/// Lower-cased file name with any `.gz` suffix removed. Every extension decision goes through
/// this, so the "is it gzipped", "is it FASTQ" and "can we digest it" predicates cannot
/// disagree — a mismatch there silently feeds compressed bytes to the line parser.
fn logical_name(path: &Path) -> String {
    let name = path.to_string_lossy().to_ascii_lowercase();
    name.strip_suffix(".gz").unwrap_or(&name).to_string()
}

/// Is this path gzipped? Case-insensitive, so `.GZ` is handled like `.gz`.
pub fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

/// Strip a trailing `.gz` (if any) and then the file extension, yielding the record/genome name.
///
/// `Path::file_stem` alone removes only the last suffix, so `GCF_000001.fna.gz` would become
/// `GCF_000001.fna` — a different identifier than the same genome uncompressed, which silently
/// breaks joins against truth tables keyed on accession.
pub fn fastx_stem(path: &Path) -> String {
    let base = if is_gzipped(path) {
        Path::new(path.file_stem().unwrap_or_default()).to_path_buf()
    } else {
        path.to_path_buf()
    };
    base.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("?")
        .to_string()
}

/// Is this a FASTQ path (after stripping any `.gz`)?
fn is_fastq(path: &Path) -> bool {
    let stem = logical_name(path);
    stem.ends_with(".fq") || stem.ends_with(".fastq")
}

/// Does this look like a **FASTA** reference we can digest (plain or gzipped)?
///
/// Deliberately excludes FASTQ: a stray read file left in a genome directory would otherwise be
/// digested as if it were a reference genome, silently adding a bogus unit to the database.
pub fn is_fasta_path(path: &Path) -> bool {
    let stem = logical_name(path);
    [".fa", ".fasta", ".fna"].iter().any(|e| stem.ends_with(e))
}

/// Does this look like a sequence file we can digest (FASTA or FASTQ, plain or gzipped)?
pub fn is_fastx_path(path: &Path) -> bool {
    let stem = logical_name(path);
    [".fa", ".fasta", ".fna", ".fq", ".fastq"]
        .iter()
        .any(|e| stem.ends_with(e))
}

/// Open a possibly-gzipped file for buffered reading.
///
/// Decompression shells out to `gzip -dc` rather than linking an inflate implementation —
/// `gzip` is present on every macOS/Linux host and this keeps the crate dependency-free.
fn open_maybe_gz(path: &Path) -> io::Result<(Box<dyn BufRead>, Option<GzipChild>)> {
    const BUF: usize = 1 << 20;
    if is_gzipped(path) {
        // Open first, so a missing/unreadable file reports ENOENT like the plain path does
        // instead of being blamed on a corrupt archive.
        let f = File::open(path)?;
        let mut child = Command::new("gzip")
            .arg("-dc")
            .stdin(Stdio::from(f))
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| {
                io::Error::other(format!(
                    "cannot decompress {}: failed to run `gzip -dc` ({e}). \
                     Install gzip or decompress the file first.",
                    path.display()
                ))
            })?;
        let out = child.stdout.take().expect("stdout was piped");
        Ok((
            Box::new(BufReader::with_capacity(BUF, out)),
            Some(GzipChild(child)),
        ))
    } else {
        let f = File::open(path)?;
        Ok((Box::new(BufReader::with_capacity(BUF, f)), None))
    }
}

/// Owns the `gzip -dc` process and reaps it on **every** exit path.
///
/// `std::process::Child` has no `Drop` that waits, so an early `?` return or a panic in the
/// caller's closure would otherwise leave a zombie behind — and `digest_genome_dir` fans out
/// over hundreds of `.fna.gz` files, so one failure per file accumulates.
struct GzipChild(Child);

impl GzipChild {
    /// Wait and turn a non-zero exit into an error (truncated/corrupt archive).
    fn finish(mut self, path: &Path) -> io::Result<()> {
        let status = self.0.wait()?;
        // Defuse the Drop guard: we have already reaped.
        std::mem::forget(self);
        if !status.success() {
            return Err(io::Error::other(format!(
                "gzip -dc failed on {} (truncated or corrupt archive?)",
                path.display()
            )));
        }
        Ok(())
    }
}

impl Drop for GzipChild {
    fn drop(&mut self) {
        // Abandoned mid-read: kill the writer so it does not block on a full pipe, then reap.
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// Reap the decompressor, if there was one.
fn reap(child: Option<GzipChild>, path: &Path) -> io::Result<()> {
    match child {
        Some(c) => c.finish(path),
        None => Ok(()),
    }
}

/// Stream every sequence in a FASTA/FASTQ file (plain or `.gz`), calling `f` on each.
///
/// Memory is bounded by the longest single record, not the file size — this is what makes
/// multi-GB samples possible. FASTA records with wrapped lines are concatenated; FASTQ takes
/// line 2 of every 4-line record.
pub fn for_each_sequence<F: FnMut(&[u8])>(path: &Path, mut f: F) -> io::Result<()> {
    let (mut reader, child) = open_maybe_gz(path)?;
    let fastq = is_fastq(path);
    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut cur: Vec<u8> = Vec::new();
    let mut lineno = 0usize;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        // Trim ASCII whitespace, not just the EOL bytes: the whole-file readers this replaced
        // used `str::trim`, and a FASTA line padded with trailing spaces must not splice them
        // into the contig (they would break `is_pure_atcg` for every window spanning them and
        // shift every downstream tag, changing marker values versus earlier versions).
        while matches!(line.last(), Some(b) if b.is_ascii_whitespace()) {
            line.pop();
        }
        let start = line
            .iter()
            .position(|b| !b.is_ascii_whitespace())
            .unwrap_or(line.len());
        let line = &line[start..];
        if fastq {
            if lineno % 4 == 1 {
                f(line);
            }
            lineno += 1;
        } else if line.first() == Some(&b'>') {
            if !cur.is_empty() {
                f(&cur);
                cur.clear();
            }
        } else {
            cur.extend_from_slice(line);
        }
    }
    if !fastq && !cur.is_empty() {
        f(&cur);
    }
    reap(child, path)
}

/// Read every sequence into memory. Use only where the whole set is genuinely needed
/// (a single reference genome); samples should go through [`sample_marker_counts_stream`].
pub fn read_fastx(path: &Path) -> io::Result<Vec<Vec<u8>>> {
    let mut seqs = Vec::new();
    for_each_sequence(path, |s| seqs.push(s.to_vec()))?;
    Ok(seqs)
}

/// Sequences per parallel digest batch. Bounds peak memory at roughly
/// `STREAM_BATCH × read_length` while keeping every worker thread fed.
pub const STREAM_BATCH: usize = 1 << 16;

/// Stream a sample and accumulate marker counts, digesting in parallel batches.
///
/// This replaces "load the whole FASTQ into a `Vec<Vec<u8>>`, then digest": peak memory is now
/// one batch (~tens of MB) instead of the entire file, so 20 GB samples run in bounded RAM.
/// Batch buffers are recycled across batches, so steady-state allocation is zero.
pub fn sample_marker_counts_stream(path: &Path, enzymes: &[&Enzyme]) -> io::Result<MarkerCounts> {
    let mut total = MarkerCounts::default();
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(STREAM_BATCH);
    let mut used = 0usize;

    for_each_sequence(path, |s| {
        // Recycle the previous batch's allocations instead of reallocating per read.
        if used < batch.len() {
            batch[used].clear();
            batch[used].extend_from_slice(s);
        } else {
            batch.push(s.to_vec());
        }
        used += 1;

        if used >= STREAM_BATCH {
            let partial = genome_marker_counts_multi_par(&batch[..used], enzymes);
            merge_counts(&mut total, partial);
            used = 0;
        }
    })?;

    if used > 0 {
        let partial = genome_marker_counts_multi_par(&batch[..used], enzymes);
        merge_counts(&mut total, partial);
    }
    Ok(total)
}

#[inline]
fn merge_counts(into: &mut MarkerCounts, from: MarkerCounts) {
    if into.is_empty() {
        *into = from;
        return;
    }
    for (m, c) in from {
        *into.entry(m).or_insert(0) += c;
    }
}

/// Digest all sequences and return the **set** of markers (for a reference genome).
pub fn genome_markers(seqs: &[Vec<u8>], enzyme: &Enzyme) -> Vec<Marker> {
    let mut set: FxHashSet<Marker> = FxHashSet::default();
    for s in seqs {
        enzyme.for_each_tag(s, |pos, len| {
            set.insert(marker_from_tag(&s[pos..pos + len]));
        });
    }
    set.into_iter().collect()
}

/// Per-genome marker copy numbers (how many times each tag occurs in the genome).
pub fn genome_marker_counts(seqs: &[Vec<u8>], enzyme: &Enzyme) -> MarkerCounts {
    genome_marker_counts_multi(seqs, &[enzyme])
}

/// Single-copy markers only (tags occurring exactly once in the genome).
///
/// StrainScan and Fast2bRAD-M (`remove_redundant`) both restrict markers to single-copy
/// loci, because multi-copy tags inflate abundance estimates and blur strain identity.
pub fn single_copy_markers(counts: &MarkerCounts) -> Vec<Marker> {
    counts
        .iter()
        .filter(|(_, &c)| c == 1)
        .map(|(&m, _)| m)
        .collect()
}

/// Digest all reads and return marker **counts** (for a sample).
pub fn sample_marker_counts(seqs: &[Vec<u8>], enzyme: &Enzyme) -> MarkerCounts {
    sample_marker_counts_multi(seqs, &[enzyme])
}

/// Per-genome marker copy numbers with a set of enzymes (pooled).
pub fn genome_marker_counts_multi(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> MarkerCounts {
    let mut counts = MarkerCounts::default();
    for s in seqs {
        count_markers_into(s, enzymes, &mut counts);
    }
    counts
}

/// Sample marker counts with a set of enzymes (pooled).
pub fn sample_marker_counts_multi(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> MarkerCounts {
    genome_marker_counts_multi(seqs, enzymes)
}

/// Parallel marker counts: digest sequence chunks across threads, then merge the maps.
/// Read digestion dominates per-sample profiling time, so this is the main speedup path.
pub fn genome_marker_counts_multi_par(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> MarkerCounts {
    let nt = crate::parallel::num_threads().min(seqs.len().max(1));
    if nt <= 1 || seqs.len() < 4096 {
        return genome_marker_counts_multi(seqs, enzymes);
    }
    let chunk = seqs.len().div_ceil(nt);
    let chunks: Vec<&[Vec<u8>]> = seqs.chunks(chunk).collect();
    let partials: Vec<MarkerCounts> =
        crate::parallel::par_map(&chunks, |c| genome_marker_counts_multi(c, enzymes));
    let mut out = MarkerCounts::default();
    for p in partials {
        merge_counts(&mut out, p);
    }
    out
}

/// Back-compat alias for [`genome_marker_counts_multi_par`].
pub fn sample_marker_counts_multi_par(seqs: &[Vec<u8>], enzymes: &[&Enzyme]) -> MarkerCounts {
    genome_marker_counts_multi_par(seqs, enzymes)
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

    /// The allocation-free `marker_from_tag` must agree with the reference definition
    /// (`hash_bytes(canonical(uppercase(tag)))`) on every input, including N-containing and
    /// palindromic tags — this is what guarantees databases stay readable across the rewrite.
    #[test]
    fn fast_path_matches_reference_canonicalization() {
        fn reference(tag: &[u8]) -> Marker {
            let mut up = tag.to_vec();
            up.make_ascii_uppercase();
            hash_bytes(&canonical(&up))
        }
        let alphabet = b"ACGTNacgtn";
        // deterministic pseudo-random tags across every enzyme tag length
        for len in [25usize, 27, 28, 32, 33] {
            for seed in 0..400u64 {
                let tag: Vec<u8> = (0..len)
                    .map(|i| {
                        let x = seed
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(i as u64 * 1442695040888963407);
                        alphabet[((x >> 33) % alphabet.len() as u64) as usize]
                    })
                    .collect();
                assert_eq!(marker_from_tag(&tag), reference(&tag), "tag {tag:?}");
            }
        }
        // explicit edge cases: palindrome, all-N, single base
        for tag in [&b"ACGT"[..], &b"NNNN"[..], &b"A"[..], &b"AT"[..], &b"GC"[..]] {
            assert_eq!(marker_from_tag(tag), reference(tag), "tag {tag:?}");
        }
    }

    #[test]
    fn digestion_is_strand_invariant() {
        // Digesting a sequence and its reverse complement must yield the SAME marker set, for
        // every enzyme. This holds because each enzyme's forward/reverse patterns are exact
        // reverse-complement pairs at its tag length, so a forward-only scan + canonical hashing
        // is strand-invariant (mirroring Fast2bRAD-M / 2bRADExtraction.pl) — no both-strand scan.
        use crate::enzymes::ALL_ENZYMES;
        let seq: Vec<u8> = (0..8000)
            .map(|i| b"ACGT"[((i * 7 + i / 3 + (i % 11)) % 4) as usize])
            .collect();
        let rc = revcomp(&seq);
        for e in ALL_ENZYMES {
            let fwd: FxHashSet<Marker> = digest_sequence(&seq, e).into_iter().collect();
            let rev: FxHashSet<Marker> = digest_sequence(&rc, e).into_iter().collect();
            assert_eq!(fwd, rev, "{} digestion must be strand-invariant", e.name);
        }
    }

    /// Digestion must be case-insensitive: soft-masked (lower-case) reference regions have to
    /// yield the same markers as the upper-case sequence, now that we scan buffers in place.
    #[test]
    fn digestion_is_case_insensitive() {
        use crate::enzymes::ALL_ENZYMES;
        let seq: Vec<u8> = (0..4000)
            .map(|i| b"ACGT"[((i * 5 + i / 7 + (i % 13)) % 4) as usize])
            .collect();
        let mut lower = seq.clone();
        lower.make_ascii_lowercase();
        for e in ALL_ENZYMES {
            assert_eq!(
                digest_sequence(&seq, e),
                digest_sequence(&lower, e),
                "{} must digest soft-masked sequence identically",
                e.name
            );
        }
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

    /// A genome must get the same identifier whether or not it was supplied gzipped.
    /// `Path::file_stem` alone leaves `X.fna.gz` named `X.fna`, which silently breaks joins
    /// against truth tables keyed on accession.
    #[test]
    fn gz_and_plain_paths_yield_the_same_stem() {
        for (p, want) in [
            ("/d/GCF_000001.fna", "GCF_000001"),
            ("/d/GCF_000001.fna.gz", "GCF_000001"),
            ("/d/GCF_000001.FNA.GZ", "GCF_000001"),
            ("/d/sample.fq.gz", "sample"),
            ("/d/a.b.c.fasta.gz", "a.b.c"),
        ] {
            assert_eq!(fastx_stem(Path::new(p)), want, "stem of {p}");
        }
    }

    /// The three extension predicates must agree on what counts as gzipped, or an uppercase
    /// `.GZ` passes the input filter and is then read as raw compressed bytes — yielding zero
    /// markers with no error at all.
    #[test]
    fn extension_predicates_agree_and_are_case_insensitive() {
        for p in ["/d/x.fna.gz", "/d/x.FNA.GZ", "/d/x.Fna.Gz"] {
            assert!(is_gzipped(Path::new(p)), "{p} must be gzipped");
            assert!(is_fasta_path(Path::new(p)), "{p} must be FASTA");
            assert!(is_fastx_path(Path::new(p)), "{p} must be FASTX");
        }
        assert!(!is_gzipped(Path::new("/d/x.fna")));
        // FASTQ is FASTX but NOT FASTA: a stray read file in a genome dir must not be
        // digested as a reference genome.
        for p in ["/d/reads.fq", "/d/reads.fastq.gz", "/d/reads.FQ.GZ"] {
            assert!(is_fastx_path(Path::new(p)), "{p} must be FASTX");
            assert!(!is_fasta_path(Path::new(p)), "{p} must NOT count as FASTA");
        }
    }

    /// FASTA lines padded with trailing whitespace must not splice spaces into the contig —
    /// that shifts every downstream tag and changes marker values versus earlier versions.
    #[test]
    fn fasta_line_whitespace_is_trimmed() {
        use std::io::Write;
        let path = std::env::temp_dir().join("s2bs_ws_test.fa");
        {
            let mut f = File::create(&path).unwrap();
            write!(f, ">a\r\nACGT  \r\n  GGGG\t\r\n").unwrap();
        }
        let seqs = read_fastx(&path).unwrap();
        assert_eq!(seqs, vec![b"ACGTGGGG".to_vec()]);
        let _ = std::fs::remove_file(&path);
    }

    /// A missing gzipped file must report "not found", not "corrupt archive".
    #[test]
    fn missing_gz_reports_not_found() {
        let err = read_fastx(Path::new("/nonexistent/dir/missing.fa.gz")).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound, "got: {err}");
    }

    /// Streaming and whole-file reading must agree, for FASTA and FASTQ, plain and gzipped.
    #[test]
    fn streaming_matches_whole_file_read() {
        use crate::enzymes::BCGI;
        use std::io::Write;

        let dir = std::env::temp_dir();
        let fa = dir.join("s2bs_stream_test.fa");
        let fq = dir.join("s2bs_stream_test.fq");

        let mut w = vec![b'A'; 32];
        w[10..13].copy_from_slice(b"CGA");
        w[19..22].copy_from_slice(b"TGC");

        // FASTA with wrapped lines
        {
            let mut f = File::create(&fa).unwrap();
            for r in 0..3 {
                writeln!(f, ">contig{r}").unwrap();
                writeln!(f, "{}", String::from_utf8(w[..16].to_vec()).unwrap()).unwrap();
                writeln!(f, "{}", String::from_utf8(w[16..].to_vec()).unwrap()).unwrap();
            }
        }
        // FASTQ
        {
            let mut f = File::create(&fq).unwrap();
            for r in 0..3 {
                writeln!(f, "@read{r}").unwrap();
                writeln!(f, "{}", String::from_utf8(w.clone()).unwrap()).unwrap();
                writeln!(f, "+").unwrap();
                writeln!(f, "{}", "I".repeat(w.len())).unwrap();
            }
        }

        for path in [&fa, &fq] {
            let seqs = read_fastx(path).unwrap();
            assert_eq!(seqs.len(), 3, "{}", path.display());
            let whole = sample_marker_counts(&seqs, &BCGI);
            let streamed = sample_marker_counts_stream(path, &[&BCGI]).unwrap();
            assert_eq!(whole, streamed, "streaming mismatch on {}", path.display());

            // gzip round-trip, if gzip is available on this host
            let gz = path.with_extension(format!(
                "{}.gz",
                path.extension().unwrap().to_str().unwrap()
            ));
            let ok = Command::new("gzip")
                .arg("-kf")
                .arg(path)
                .status()
                .map(|s| s.success())
                .unwrap_or(false);
            if ok && gz.exists() {
                let from_gz = sample_marker_counts_stream(&gz, &[&BCGI]).unwrap();
                assert_eq!(whole, from_gz, "gzip mismatch on {}", gz.display());
                let _ = std::fs::remove_file(&gz);
            }
            let _ = std::fs::remove_file(path);
        }
    }
}
