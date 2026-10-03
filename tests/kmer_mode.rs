//! End-to-end k-mer-mode test: two synthetic genomes (shared core + private regions),
//! k=15, scale=1 → cluster DB (with the `kmer<K>s<S>` header token) → save/load → profile
//! a ~70/30 read mixture → both clusters recovered with plausible abundances; plus a
//! `--layer1 cst --layer2 enet` smoke run on the same DB.

use strain2bscan::cst::{SpeciesCst, DEFAULT_SIMILARITY};
use strain2bscan::db::StrainDb;
use strain2bscan::identify::{profile, Layer1, Layer2, Params};
use strain2bscan::markers::{
    genome_kmer_counts, kmer_db_token, parse_kmer_db_token, single_copy_markers, sketch_threshold,
    Marker, MarkerCounts,
};

const K: usize = 15;
const SCALE: u64 = 1;

struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn dna(&mut self, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Two ~6.5 kb genomes: a 4 kb shared core plus 2.5 kb of private sequence each, so their
/// k-mer Jaccard is ~0.4 and they must form two clusters.
fn synth_genomes() -> (Vec<u8>, Vec<u8>) {
    let mut rng = XorShift(0x9e3779b97f4a7c15);
    let core = rng.dna(4000);
    let mut a = core.clone();
    a.extend(rng.dna(2500));
    let mut b = core;
    b.extend(rng.dna(2500));
    (a, b)
}

/// Exact-weight read mixture: `n_reads[i]` random 150-bp windows of genome `i`.
fn synth_reads(genomes: &[&Vec<u8>], n_reads: &[usize], read_len: usize) -> Vec<Vec<u8>> {
    let mut rng = XorShift(0xdeadbeefcafe);
    let mut reads = Vec::new();
    for (&g, &n) in genomes.iter().zip(n_reads) {
        for _ in 0..n {
            let start = rng.below(g.len() - read_len + 1);
            reads.push(g[start..start + read_len].to_vec());
        }
    }
    reads
}

/// What `strain2bscan cluster --marker-source kmer --kmer-size 15 --sketch-scale 1` does,
/// in memory: digest → single-copy markers → CST → cluster DB + header token + tree.
fn build_kmer_cluster_db(genomes: &[Vec<u8>]) -> StrainDb {
    let threshold = sketch_threshold(SCALE);
    let recs: Vec<(String, Vec<Marker>, Vec<Marker>)> = genomes
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let counts = genome_kmer_counts(std::slice::from_ref(g), K, threshold);
            let full: Vec<Marker> = counts.keys().copied().collect();
            (format!("g{i}"), single_copy_markers(&counts), full)
        })
        .collect();
    let cst = SpeciesCst::build(recs, DEFAULT_SIMILARITY, false);
    assert_eq!(
        cst.n_clusters(),
        2,
        "core+private genomes must split into 2 clusters"
    );
    let mut db = cst.cluster_db();
    db.enzymes = vec![kmer_db_token(K, SCALE)];
    db.tree = Some(cst.build_tree());
    db
}

fn mixture_counts() -> MarkerCounts {
    let (a, b) = synth_genomes();
    let reads = synth_reads(&[&a, &b], &[700, 300], 150);
    genome_kmer_counts(&reads, K, sketch_threshold(SCALE))
}

#[test]
fn kmer_mode_end_to_end_recovers_70_30_mixture() {
    let (a, b) = synth_genomes();
    let db = build_kmer_cluster_db(&[a, b]);

    // Save → load round trip: the kmer<K>s<S> token, the tree, and the per-strain marker
    // counts (4th header field) must all survive, and the token must parse back.
    let path = std::env::temp_dir().join(format!("s2bs_kmer_e2e_{}.tsv", std::process::id()));
    db.save(&path).unwrap();
    let text = String::from_utf8(std::fs::read(&path).unwrap()).unwrap();
    let header = text.lines().next().unwrap();
    assert_eq!(
        header.split('\t').count(),
        4,
        "header must keep the counts field"
    );
    assert_eq!(header.split('\t').nth(1), Some("2"));
    assert_eq!(header.split('\t').nth(2), Some("kmer15s1"));
    let db = StrainDb::load(&path).unwrap();
    let _ = std::fs::remove_file(&path);
    assert_eq!(parse_kmer_db_token(&db.enzymes[0]), Some((K, SCALE)));
    assert!(db.tree.is_some(), "the CST must survive the round trip");

    // Profile the 70/30 mixture with default (unique-marker) parameters.
    let counts = mixture_counts();
    let calls = profile(&db, &counts, &Params::default());
    assert_eq!(calls.len(), 2, "both clusters must be called: {calls:?}");
    // Cluster order follows genome order (g0 → C0), so C0 is the 70% genome.
    let c0 = calls.iter().find(|c| c.name == "C0").unwrap();
    let c1 = calls.iter().find(|c| c.name == "C1").unwrap();
    assert!(
        (0.55..=0.85).contains(&c0.rel_abundance),
        "C0 abundance {:.3} far from 0.70",
        c0.rel_abundance
    );
    assert!(
        (0.15..=0.45).contains(&c1.rel_abundance),
        "C1 abundance {:.3} far from 0.30",
        c1.rel_abundance
    );
    let sum: f64 = calls.iter().map(|c| c.rel_abundance).sum();
    assert!(
        (sum - 1.0).abs() < 1e-6,
        "abundances must sum to 1, got {sum}"
    );
    // Both are fully covered: 700 reads of 150 bp over ~6.5 kb is ~15x per k-mer.
    assert!(c0.coverage > 0.9 && c1.coverage > 0.9);
}

#[test]
fn kmer_mode_cst_enet_smoke() {
    let (a, b) = synth_genomes();
    let db = build_kmer_cluster_db(&[a, b]);
    assert!(db.tree.is_some());
    let counts = mixture_counts();
    let params = Params {
        layer1: Layer1::Cst,
        layer2: Layer2::Enet,
        ..Params::default()
    };
    // Must not panic, and the output must be well-formed (parses as a prediction set).
    let calls = profile(&db, &counts, &params);
    assert!(
        !calls.is_empty(),
        "cst/enet on the k-mer DB must call something"
    );
    for c in &calls {
        assert!(c.rel_abundance.is_finite() && c.rel_abundance > 0.0);
        assert!(c.depth.is_finite() && c.depth >= 0.0);
        assert!((0.0..=1.0).contains(&c.coverage));
    }
    let sum: f64 = calls.iter().map(|c| c.rel_abundance).sum();
    assert!(
        (sum - 1.0).abs() < 1e-6,
        "abundances must sum to 1, got {sum}"
    );
}
