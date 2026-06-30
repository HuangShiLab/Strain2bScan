//! Type-IIB restriction enzyme digestion, ported from Fast2bRAD-M `src/enzymes.rs`.
//!
//! A type-IIB enzyme cuts on both sides of its recognition site, releasing a short
//! fixed-length fragment (the "2bRAD tag", 32–38 bp). We model each enzyme as a set of
//! `Pattern`s; a pattern is a list of `Anchor`s (a motif required at a fixed offset
//! inside a `tag_length` window). Scanning every offset of a sequence and testing the
//! anchors reproduces the enzyme's digestion sites. Forward + reverse patterns let us scan
//! only the forward strand and canonicalize.
//!
//! All 16 enzymes from Fast2bRAD-M's table are ported. Use cases:
//!   * **BcgI 2bRAD data** → digest with `BcgI` only.
//!   * **Conventional metagenome (150 bp / long reads)** → digest with `all` 16 enzymes to
//!     enrich strain-specific markers (~16× more tag loci). The genome DB must be built with
//!     the same enzyme set as the sample.
//!
//! NOTE: HaeIV and Hin4I have highly degenerate sites; Fast2bRAD's table encodes only their
//! 2 most-conserved bases, so they over-match (many tags). They are included in `all` for
//! fidelity with Fast2bRAD; drop them if they add noise on your data.

#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub offset: usize,
    pub motif: &'static [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct Pattern {
    pub anchors: &'static [Anchor],
}

impl Pattern {
    pub fn matches(&self, window: &[u8]) -> bool {
        self.anchors.iter().all(|a| {
            let end = a.offset + a.motif.len();
            end <= window.len() && &window[a.offset..end] == a.motif
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
    /// Return all distinct `(offset, tag_length)` digestion sites in `sequence`.
    pub fn find_all_tags(&self, sequence: &[u8]) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        if sequence.len() < self.tag_length {
            return out;
        }
        let last = sequence.len() - self.tag_length;
        for offset in 0..=last {
            let window = &sequence[offset..offset + self.tag_length];
            if self.patterns.iter().any(|p| p.matches(window)) && is_pure_atcg(window) {
                out.push((offset, self.tag_length));
            }
        }
        out
    }
}

// ===== Enzyme table — all 16, ported verbatim from Fast2bRAD-M =================
// (concat/pear params from Fast2bRAD are Type-4-only and omitted here.)

macro_rules! anchors {
    ($($off:expr => $motif:expr),+ $(,)?) => {
        &[$(Anchor { offset: $off, motif: $motif }),+]
    };
}

// 1. CspCI (36)
const CSPCI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(11 => b"CAA", 19 => b"GTGG"),
    },
    Pattern {
        anchors: anchors!(10 => b"CCAC", 19 => b"TTG"),
    },
];
pub const CSPCI: Enzyme = Enzyme {
    name: "CspCI",
    id: 1,
    tag_length: 36,
    patterns: &CSPCI_P,
};

// 2. AloI (37)
const ALOI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(7 => b"GAAC", 17 => b"TCC"),
    },
    Pattern {
        anchors: anchors!(7 => b"GGA", 16 => b"GTTC"),
    },
];
pub const ALOI: Enzyme = Enzyme {
    name: "AloI",
    id: 2,
    tag_length: 37,
    patterns: &ALOI_P,
};

// 3. BsaXI (32)
const BSAXI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(9 => b"AC", 16 => b"CTCC"),
    },
    Pattern {
        anchors: anchors!(7 => b"GGAG", 16 => b"GT"),
    },
];
pub const BSAXI: Enzyme = Enzyme {
    name: "BsaXI",
    id: 3,
    tag_length: 32,
    patterns: &BSAXI_P,
};

// 4. BaeI (36)
const BAEI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(10 => b"AC", 16 => b"GTA"),
    },
    Pattern {
        anchors: anchors!(7 => b"G", 9 => b"TAC"),
    },
];
pub const BAEI: Enzyme = Enzyme {
    name: "BaeI",
    id: 4,
    tag_length: 36,
    patterns: &BAEI_P,
};

// 5. BcgI (32) — the canonical 2bRAD enzyme
const BCGI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(10 => b"CGA", 19 => b"TGC"),
    },
    Pattern {
        anchors: anchors!(10 => b"GCA", 19 => b"TCG"),
    },
];
pub const BCGI: Enzyme = Enzyme {
    name: "BcgI",
    id: 5,
    tag_length: 32,
    patterns: &BCGI_P,
};

// 6. CjeI (37)
const CJEI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(8 => b"CCA", 17 => b"GT"),
    },
    Pattern {
        anchors: anchors!(9 => b"AC", 17 => b"TGG"),
    },
];
pub const CJEI: Enzyme = Enzyme {
    name: "CjeI",
    id: 6,
    tag_length: 37,
    patterns: &CJEI_P,
};

// 7. PpiI (35)
const PPII_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(7 => b"GAAC", 17 => b"CTC"),
    },
    Pattern {
        anchors: anchors!(8 => b"GAG", 16 => b"GTTC"),
    },
];
pub const PPII: Enzyme = Enzyme {
    name: "PpiI",
    id: 7,
    tag_length: 35,
    patterns: &PPII_P,
};

// 8. PsrI (35)
const PSRI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(7 => b"GAAC", 17 => b"TAC"),
    },
    Pattern {
        anchors: anchors!(7 => b"GTA", 16 => b"GTTC"),
    },
];
pub const PSRI: Enzyme = Enzyme {
    name: "PsrI",
    id: 8,
    tag_length: 35,
    patterns: &PSRI_P,
};

// 9. BplI (35, palindrome)
const BPLI_P: [Pattern; 1] = [Pattern {
    anchors: anchors!(8 => b"GAG", 16 => b"CTC"),
}];
pub const BPLI: Enzyme = Enzyme {
    name: "BplI",
    id: 9,
    tag_length: 35,
    patterns: &BPLI_P,
};

// 10. FalI (36, palindrome)
const FALI_P: [Pattern; 1] = [Pattern {
    anchors: anchors!(8 => b"AAG", 16 => b"CTT"),
}];
pub const FALI: Enzyme = Enzyme {
    name: "FalI",
    id: 10,
    tag_length: 36,
    patterns: &FALI_P,
};

// 11. Bsp24I (36)
const BSP24I_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(8 => b"GAC", 17 => b"TGG"),
    },
    Pattern {
        anchors: anchors!(7 => b"CCA", 16 => b"GTC"),
    },
];
pub const BSP24I: Enzyme = Enzyme {
    name: "Bsp24I",
    id: 11,
    tag_length: 36,
    patterns: &BSP24I_P,
};

// 12. HaeIV (37) — degenerate (single conserved dinucleotide per pattern)
const HAEIV_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(7 => b"GA"),
    },
    Pattern {
        anchors: anchors!(9 => b"GA"),
    },
];
pub const HAEIV: Enzyme = Enzyme {
    name: "HaeIV",
    id: 12,
    tag_length: 37,
    patterns: &HAEIV_P,
};

// 13. CjePI (38)
const CJEPI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(7 => b"CCA", 17 => b"TC"),
    },
    Pattern {
        anchors: anchors!(8 => b"GA", 17 => b"TGG"),
    },
];
pub const CJEPI: Enzyme = Enzyme {
    name: "CjePI",
    id: 13,
    tag_length: 38,
    patterns: &CJEPI_P,
};

// 14. Hin4I (35) — degenerate
const HIN4I_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(8 => b"GA"),
    },
    Pattern {
        anchors: anchors!(8 => b"GA"),
    },
];
pub const HIN4I: Enzyme = Enzyme {
    name: "Hin4I",
    id: 14,
    tag_length: 35,
    patterns: &HIN4I_P,
};

// 15. AlfI (33, palindrome)
const ALFI_P: [Pattern; 1] = [Pattern {
    anchors: anchors!(10 => b"GCA", 19 => b"TGC"),
}];
pub const ALFI: Enzyme = Enzyme {
    name: "AlfI",
    id: 15,
    tag_length: 33,
    patterns: &ALFI_P,
};

// 16. BslFI (33)
const BSLFI_P: [Pattern; 2] = [
    Pattern {
        anchors: anchors!(6 => b"GGGAC"),
    },
    Pattern {
        anchors: anchors!(14 => b"GTCCC"),
    },
];
pub const BSLFI: Enzyme = Enzyme {
    name: "BslFI",
    id: 16,
    tag_length: 33,
    patterns: &BSLFI_P,
};

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
    fn finds_a_crafted_bcgi_site() {
        let mut w = vec![b'A'; 32];
        w[10..13].copy_from_slice(b"CGA");
        w[19..22].copy_from_slice(b"TGC");
        let mut seq = vec![b'T'; 5];
        seq.extend_from_slice(&w);
        seq.extend_from_slice(&[b'T'; 5]);
        let hits = BCGI.find_all_tags(&seq);
        assert!(
            hits.iter().any(|&(off, len)| off == 5 && len == 32),
            "hits: {hits:?}"
        );
    }

    #[test]
    fn rejects_window_with_n() {
        let mut w = vec![b'A'; 32];
        w[10..13].copy_from_slice(b"CGA");
        w[19..22].copy_from_slice(b"TGC");
        w[0] = b'N';
        assert!(BCGI.find_all_tags(&w).is_empty());
    }

    #[test]
    fn parse_set_all_and_list() {
        assert_eq!(parse_enzyme_set("all").unwrap().len(), 16);
        let set = parse_enzyme_set("BcgI,1").unwrap();
        assert_eq!(set.len(), 2);
        assert_eq!(set[0].name, "BcgI");
        assert_eq!(set[1].name, "CspCI");
        assert!(parse_enzyme_set("nope").is_none());
    }
}
