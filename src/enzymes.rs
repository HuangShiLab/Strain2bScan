//! Type-IIB restriction enzyme digestion, ported from **Fast2bRAD-M** (`src/enzymes.rs`,
//! HuangShiLab/Fast2bRAD-M), whose tag lengths and recognition patterns are derived directly
//! from the `@site` regexes in the original `2bRADExtraction.pl` (shihuang047/2bRAD-M). Extracted
//! tags therefore match the reference Perl / Fast2bRAD-M implementation exactly for every enzyme,
//! which is required for interoperability (Fast2bRAD-M species layer → Strain2bScan strain layer)
//! and for real BcgI 2bRAD experimental data.
//!
//! A type-IIB enzyme cuts on both sides of its recognition site, releasing a short fixed-length
//! fragment (the "2bRAD tag"). Each enzyme is a set of `Pattern`s; a pattern is a list of exact
//! `Anchor`s (a literal motif required at a fixed offset in a `tag_length` window) plus, for the
//! three IUPAC-degenerate enzymes (BaeI, HaeIV, Hin4I), `ClassAnchor`s (a single position that may
//! be any of a small set of bases, e.g. `[CT]`). Non-degenerate positions are unconstrained ACGT
//! (enforced by `is_pure_atcg`), so they need no anchor.
//!
//! **Strand handling.** We scan only the forward strand, but each enzyme carries a *forward* and a
//! *reverse* pattern that are exact reverse-complement pairs at the (correct) tag length. Together
//! they find sites in either orientation in one pass, and `markers::marker_from_tag` canonicalizes
//! each tag (min of tag / its reverse complement). This makes digestion **strand-invariant**:
//! `digest(seq) == digest(revcomp(seq))` for every enzyme (see the regression test in `markers`),
//! so genome databases and reads (2bRAD tags or shotgun fragments, in either orientation) live in
//! one marker space — with no need to scan both strands.
//!
//! All 16 enzymes from the Fast2bRAD-M table are ported. HaeIV and Hin4I are highly degenerate and
//! excluded from the recommended multi-enzyme set (noisy), but supported for fidelity.

#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub offset: usize,
    /// literal motif required at `offset`
    pub motif: &'static [u8],
}

/// A single degenerate position: `window[offset]` must be one of `allowed` (IUPAC class).
#[derive(Debug, Clone, Copy)]
pub struct ClassAnchor {
    pub offset: usize,
    pub allowed: &'static [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct Pattern {
    pub anchors: &'static [Anchor],
    pub classes: &'static [ClassAnchor],
}

impl Pattern {
    /// Match `window` against this pattern, **case-insensitively**.
    ///
    /// Case insensitivity here is what lets the digest hot path scan the caller's sequence
    /// buffer directly: without it every read/contig had to be copied and upper-cased before
    /// scanning (one heap allocation + one full pass per sequence). Soft-masked reference
    /// genomes (lower-case repeat regions) are also handled correctly as a result.
    #[inline]
    pub fn matches(&self, window: &[u8]) -> bool {
        self.anchors.iter().all(|a| {
            let end = a.offset + a.motif.len();
            end <= window.len() && window[a.offset..end].eq_ignore_ascii_case(a.motif)
        }) && self.classes.iter().all(|c| {
            c.offset < window.len() && c.allowed.contains(&window[c.offset].to_ascii_uppercase())
        })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Enzyme {
    pub name: &'static str,
    pub id: u8,
    pub tag_length: usize,
    pub patterns: &'static [Pattern],
}

/// O(1) ATCG membership test (upper/lower case); everything else (incl. N) is false.
const ATCG_TABLE: [bool; 256] = {
    let mut t = [false; 256];
    t[b'A' as usize] = true;
    t[b'a' as usize] = true;
    t[b'T' as usize] = true;
    t[b't' as usize] = true;
    t[b'C' as usize] = true;
    t[b'c' as usize] = true;
    t[b'G' as usize] = true;
    t[b'g' as usize] = true;
    t
};

#[inline]
fn is_pure_atcg(window: &[u8]) -> bool {
    window.iter().all(|&b| ATCG_TABLE[b as usize])
}

impl Enzyme {
    /// Call `f(offset, tag_length)` for every digestion site in `sequence` (forward-strand
    /// scan; the forward+reverse patterns catch both orientations).
    ///
    /// This is the allocation-free form used by the digest hot path — the previous
    /// `find_all_tags` built a `Vec` per enzyme per sequence, i.e. up to 16 allocations for
    /// every read in a multi-enzyme digest.
    #[inline]
    pub fn for_each_tag<F: FnMut(usize, usize)>(&self, sequence: &[u8], mut f: F) {
        if sequence.len() < self.tag_length {
            return;
        }
        let last = sequence.len() - self.tag_length;
        for offset in 0..=last {
            let window = &sequence[offset..offset + self.tag_length];
            if self.patterns.iter().any(|p| p.matches(window)) && is_pure_atcg(window) {
                f(offset, self.tag_length);
            }
        }
    }

    /// Collecting form of [`Enzyme::for_each_tag`], kept for tests and external callers.
    pub fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        self.for_each_tag(sequence, |pos, len| out.push((pos, len)));
        out
    }
}

// ===== Enzyme table — all 16, tag lengths + patterns from Fast2bRAD-M / 2bRADExtraction.pl =====
//
// Each enzyme's forward and reverse patterns are exact reverse-complement pairs at its tag length
// (self-RC single pattern for the palindromic enzymes BplI/FalI/AlfI), so a forward-only scan +
// canonical hashing is strand-invariant.

macro_rules! anchors {
    ($($off:expr => $motif:expr),+ $(,)?) => {
        &[$(Anchor { offset: $off, motif: $motif }),+]
    };
}
macro_rules! classes {
    ($($off:expr => $set:expr),+ $(,)?) => {
        &[$(ClassAnchor { offset: $off, allowed: $set }),+]
    };
}
/// Exact (no degenerate positions) pattern.
const fn pat(anchors: &'static [Anchor]) -> Pattern {
    Pattern { anchors, classes: &[] }
}

// 1. CspCI (33)
const CSPCI_P: [Pattern; 2] = [
    pat(anchors!(11 => b"CAA", 19 => b"GTGG")),
    pat(anchors!(10 => b"CCAC", 19 => b"TTG")),
];
pub const CSPCI: Enzyme = Enzyme { name: "CspCI", id: 1, tag_length: 33, patterns: &CSPCI_P };

// 2. AloI (27)
const ALOI_P: [Pattern; 2] = [
    pat(anchors!(7 => b"GAAC", 17 => b"TCC")),
    pat(anchors!(7 => b"GGA", 16 => b"GTTC")),
];
pub const ALOI: Enzyme = Enzyme { name: "AloI", id: 2, tag_length: 27, patterns: &ALOI_P };

// 3. BsaXI (27)
const BSAXI_P: [Pattern; 2] = [
    pat(anchors!(9 => b"AC", 16 => b"CTCC")),
    pat(anchors!(7 => b"GGAG", 16 => b"GT")),
];
pub const BSAXI: Enzyme = Enzyme { name: "BsaXI", id: 3, tag_length: 27, patterns: &BSAXI_P };

// 4. BaeI (28) — degenerate: forward [CT]@19, reverse [AG]@8
const BAEI_P: [Pattern; 2] = [
    Pattern { anchors: anchors!(10 => b"AC", 16 => b"GTA", 20 => b"C"), classes: classes!(19 => b"CT") },
    Pattern { anchors: anchors!(7 => b"G", 9 => b"TAC", 16 => b"GT"), classes: classes!(8 => b"AG") },
];
pub const BAEI: Enzyme = Enzyme { name: "BaeI", id: 4, tag_length: 28, patterns: &BAEI_P };

// 5. BcgI (32) — the canonical 2bRAD enzyme
const BCGI_P: [Pattern; 2] = [
    pat(anchors!(10 => b"CGA", 19 => b"TGC")),
    pat(anchors!(10 => b"GCA", 19 => b"TCG")),
];
pub const BCGI: Enzyme = Enzyme { name: "BcgI", id: 5, tag_length: 32, patterns: &BCGI_P };

// 6. CjeI (28)
const CJEI_P: [Pattern; 2] = [
    pat(anchors!(8 => b"CCA", 17 => b"GT")),
    pat(anchors!(9 => b"AC", 17 => b"TGG")),
];
pub const CJEI: Enzyme = Enzyme { name: "CjeI", id: 6, tag_length: 28, patterns: &CJEI_P };

// 7. PpiI (27)
const PPII_P: [Pattern; 2] = [
    pat(anchors!(7 => b"GAAC", 16 => b"CTC")),
    pat(anchors!(8 => b"GAG", 16 => b"GTTC")),
];
pub const PPII: Enzyme = Enzyme { name: "PpiI", id: 7, tag_length: 27, patterns: &PPII_P };

// 8. PsrI (27)
const PSRI_P: [Pattern; 2] = [
    pat(anchors!(7 => b"GAAC", 17 => b"TAC")),
    pat(anchors!(7 => b"GTA", 16 => b"GTTC")),
];
pub const PSRI: Enzyme = Enzyme { name: "PsrI", id: 8, tag_length: 27, patterns: &PSRI_P };

// 9. BplI (27, palindrome — self-RC single pattern)
const BPLI_P: [Pattern; 1] = [pat(anchors!(8 => b"GAG", 16 => b"CTC"))];
pub const BPLI: Enzyme = Enzyme { name: "BplI", id: 9, tag_length: 27, patterns: &BPLI_P };

// 10. FalI (27, palindrome)
const FALI_P: [Pattern; 1] = [pat(anchors!(8 => b"AAG", 16 => b"CTT"))];
pub const FALI: Enzyme = Enzyme { name: "FalI", id: 10, tag_length: 27, patterns: &FALI_P };

// 11. Bsp24I (27)
const BSP24I_P: [Pattern; 2] = [
    pat(anchors!(8 => b"GAC", 17 => b"TGG")),
    pat(anchors!(7 => b"CCA", 16 => b"GTC")),
];
pub const BSP24I: Enzyme = Enzyme { name: "Bsp24I", id: 11, tag_length: 27, patterns: &BSP24I_P };

// 12. HaeIV (27) — degenerate: forward [CT]@9,[AG]@15 ; reverse [CT]@11,[AG]@17
const HAEIV_P: [Pattern; 2] = [
    Pattern { anchors: anchors!(7 => b"GA", 16 => b"TC"), classes: classes!(9 => b"CT", 15 => b"AG") },
    Pattern { anchors: anchors!(9 => b"GA", 18 => b"TC"), classes: classes!(11 => b"CT", 17 => b"AG") },
];
pub const HAEIV: Enzyme = Enzyme { name: "HaeIV", id: 12, tag_length: 27, patterns: &HAEIV_P };

// 13. CjePI (27)
const CJEPI_P: [Pattern; 2] = [
    pat(anchors!(7 => b"CCA", 17 => b"TC")),
    pat(anchors!(8 => b"GA", 17 => b"TGG")),
];
pub const CJEPI: Enzyme = Enzyme { name: "CjePI", id: 13, tag_length: 27, patterns: &CJEPI_P };

// 14. Hin4I (27) — degenerate: forward [CT]@10,[GAC]@16 ; reverse [CTG]@10,[AG]@16
const HIN4I_P: [Pattern; 2] = [
    Pattern { anchors: anchors!(8 => b"GA", 17 => b"TC"), classes: classes!(10 => b"CT", 16 => b"GAC") },
    Pattern { anchors: anchors!(8 => b"GA", 17 => b"TC"), classes: classes!(10 => b"CTG", 16 => b"AG") },
];
pub const HIN4I: Enzyme = Enzyme { name: "Hin4I", id: 14, tag_length: 27, patterns: &HIN4I_P };

// 15. AlfI (32, palindrome)
const ALFI_P: [Pattern; 1] = [pat(anchors!(10 => b"GCA", 19 => b"TGC"))];
pub const ALFI: Enzyme = Enzyme { name: "AlfI", id: 15, tag_length: 32, patterns: &ALFI_P };

// 16. BslFI (25)
const BSLFI_P: [Pattern; 2] = [
    pat(anchors!(6 => b"GGGAC")),
    pat(anchors!(14 => b"GTCCC")),
];
pub const BSLFI: Enzyme = Enzyme { name: "BslFI", id: 16, tag_length: 25, patterns: &BSLFI_P };

/// All 16 enzymes, in id order.
pub static ALL_ENZYMES: &[&Enzyme] = &[
    &CSPCI, &ALOI, &BSAXI, &BAEI, &BCGI, &CJEI, &PPII, &PSRI, &BPLI, &FALI, &BSP24I, &HAEIV,
    &CJEPI, &HIN4I, &ALFI, &BSLFI,
];

pub fn enzyme_by_name(name: &str) -> Option<&'static Enzyme> {
    ALL_ENZYMES
        .iter()
        .copied()
        .find(|e| e.name.eq_ignore_ascii_case(name))
}

pub fn enzyme_by_id(id: u8) -> Option<&'static Enzyme> {
    ALL_ENZYMES.iter().copied().find(|e| e.id == id)
}

/// Accept either an enzyme name ("BcgI") or its numeric id ("5").
pub fn parse_enzyme(site: &str) -> Option<&'static Enzyme> {
    enzyme_by_name(site).or_else(|| site.parse::<u8>().ok().and_then(enzyme_by_id))
}

/// Parse an enzyme *set* spec: `all` (all 16), or a comma list of names/ids ("BcgI,CspCI").
pub fn parse_enzyme_set(spec: &str) -> Option<Vec<&'static Enzyme>> {
    if spec.eq_ignore_ascii_case("all") {
        return Some(ALL_ENZYMES.to_vec());
    }
    let mut out = Vec::new();
    for tok in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        out.push(parse_enzyme(tok)?);
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_sixteen_present_with_unique_ids() {
        assert_eq!(ALL_ENZYMES.len(), 16);
        let mut ids: Vec<u8> = ALL_ENZYMES.iter().map(|e| e.id).collect();
        ids.sort();
        assert_eq!(ids, (1..=16).collect::<Vec<u8>>());
    }

    #[test]
    fn tag_lengths_match_fast2brad_m() {
        // From the authoritative Fast2bRAD-M / 2bRADExtraction.pl table.
        let expected: [(&str, usize); 16] = [
            ("CspCI", 33), ("AloI", 27), ("BsaXI", 27), ("BaeI", 28), ("BcgI", 32), ("CjeI", 28),
            ("PpiI", 27), ("PsrI", 27), ("BplI", 27), ("FalI", 27), ("Bsp24I", 27), ("HaeIV", 27),
            ("CjePI", 27), ("Hin4I", 27), ("AlfI", 32), ("BslFI", 25),
        ];
        for (name, len) in expected {
            let e = enzyme_by_name(name).unwrap();
            assert_eq!(e.tag_length, len, "{name} tag_length");
        }
    }

    #[test]
    fn finds_a_crafted_bcgi_site() {
        let mut w = vec![b'A'; 32];
        w[10..13].copy_from_slice(b"CGA");
        w[19..22].copy_from_slice(b"TGC");
        let mut seq = vec![b'T'; 5];
        seq.extend_from_slice(&w);
        seq.extend_from_slice(&[b'T'; 5]);
        let hits = BCGI.find_all_tags(&seq);
        assert!(hits.iter().any(|&(pos, len)| pos == 5 && len == 32));
    }

    #[test]
    fn forward_reverse_patterns_are_rc_pairs() {
        // Each enzyme's patterns must be exact reverse-complement pairs at its tag length, which
        // is what makes forward-only scan + canonical hashing strand-invariant. We check this by
        // building a window that matches the forward pattern and confirming its reverse complement
        // matches (some) pattern of the same enzyme.
        fn rc(seq: &[u8]) -> Vec<u8> {
            seq.iter().rev().map(|&b| match b {
                b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', x => x,
            }).collect()
        }
        for e in ALL_ENZYMES {
            for p in e.patterns {
                // materialize a concrete window satisfying pattern p (fill gaps with 'A')
                let mut w = vec![b'A'; e.tag_length];
                for a in p.anchors {
                    w[a.offset..a.offset + a.motif.len()].copy_from_slice(a.motif);
                }
                for c in p.classes {
                    w[c.offset] = c.allowed[0];
                }
                assert!(p.matches(&w), "{} self-match", e.name);
                let r = rc(&w);
                assert!(
                    e.patterns.iter().any(|q| q.matches(&r)),
                    "{}: reverse complement of a forward-pattern window is not matched by any pattern \
                     (patterns are not RC pairs at tag_length {})",
                    e.name, e.tag_length
                );
            }
        }
    }
}
