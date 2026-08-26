//! End-to-end batch-mode test: two synthetic two-genome species DBs (BcgI enzyme mode),
//! two samples — one paired-end, one single-end — profiled through the real
//! `strain2bscan batch` binary, and every merged row checked field-by-field against a
//! separate `multi-profile --out` run on the same sample. The paired sample's reference run
//! digests `cat R1 R2`, pinning the R1+R2 count-merging equivalence.
//!
//! Same fixture style as `kmer_mode.rs` (synthetic genomes + exact-weight read windows),
//! but the commands run as subprocesses via `CARGO_BIN_EXE_*` because `batch`/`multi-profile`
//! live in the binary, not the library.

use std::path::{Path, PathBuf};
use std::process::Command;

use strain2bscan::db::StrainDb;
use strain2bscan::enzymes::parse_enzyme_set;
use strain2bscan::markers::{genome_marker_counts_multi, single_copy_markers};

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
        (0..len).map(|_| b"ACGT"[(self.next() % 4) as usize]).collect()
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Two ~200 kb genomes of one species: a 120 kb shared core plus 80 kb of private sequence
/// each, so both carry strain-unique BcgI tags. BcgI's 6 anchored bases in a 32 bp window
/// give ~1 tag per 2 kb, i.e. ~100 tags per genome — enough for the gates at
/// `--min-support 1`.
fn synth_species(seed: u64) -> (Vec<u8>, Vec<u8>) {
    let mut rng = XorShift(seed);
    let core = rng.dna(120_000);
    let mut a = core.clone();
    a.extend(rng.dna(80_000));
    let mut b = core;
    b.extend(rng.dna(80_000));
    (a, b)
}

/// What `strain2bscan build --enzyme BcgI` does, in memory: digest → single-copy markers →
/// flat DB (no clustering — nothing here needs a tree).
fn build_species_db(genomes: &[Vec<u8>], prefix: &str) -> StrainDb {
    let set = parse_enzyme_set("BcgI").unwrap();
    let strains = genomes
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let counts = genome_marker_counts_multi(std::slice::from_ref(g), &set);
            (format!("{prefix}{i}"), single_copy_markers(&counts))
        })
        .collect();
    let mut db = StrainDb::build(strains);
    db.enzymes = vec!["BcgI".to_string()];
    db
}

/// Exact-weight read mixture: `n_reads[i]` random 150-bp windows of genome `i`.
fn synth_reads(genomes: &[&Vec<u8>], n_reads: &[usize], seed: u64, read_len: usize) -> Vec<Vec<u8>> {
    let mut rng = XorShift(seed);
    let mut reads = Vec::new();
    for (&g, &n) in genomes.iter().zip(n_reads) {
        for _ in 0..n {
            let start = rng.below(g.len() - read_len + 1);
            reads.push(g[start..start + read_len].to_vec());
        }
    }
    reads
}

fn write_fastq(path: &Path, reads: &[Vec<u8>]) {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for (i, r) in reads.iter().enumerate() {
        writeln!(w, "@r{i}").unwrap();
        w.write_all(r).unwrap();
        writeln!(w, "\n+").unwrap();
        writeln!(w, "{}", "I".repeat(r.len())).unwrap();
    }
}

fn run(args: &[&str]) -> std::process::Output {
    let out = Command::new(env!("CARGO_BIN_EXE_strain2bscan"))
        .args(args)
        .output()
        .unwrap();
    out
}

fn run_ok(args: &[&str]) -> std::process::Output {
    let out = run(args);
    assert!(
        out.status.success(),
        "strain2bscan {:?} failed: {}",
        &args[..1],
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// Identification flags shared by every run, with the gates lowered to fit the tiny
/// synthetic panels (the defaults are calibrated on real ~200-marker species panels).
/// Passing them here also exercises batch's parameter passthrough.
const GATES: &[&str] = &[
    "--enzyme",
    "BcgI",
    "--min-species-markers",
    "1",
    "--min-species-detect",
    "1",
    "--min-support",
    "1",
    "--min-coverage",
    "0",
    "--min-abundance",
    "0",
    "--min-consistency",
    "0",
    "--fixed-gate",
];

/// Append the shared gate flags to a fixed argument list.
fn with_gates(mut args: Vec<String>) -> Vec<String> {
    args.extend(GATES.iter().map(|s| s.to_string()));
    args
}

fn run_strings(args: &[String]) -> std::process::Output {
    run_ok(&args.iter().map(String::as_str).collect::<Vec<_>>())
}

/// Data lines (header dropped) of a prediction TSV.
fn body_lines(path: &Path) -> Vec<String> {
    String::from_utf8(std::fs::read(path).unwrap())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

struct Fixture {
    dir: PathBuf,
}

impl Fixture {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("s2bs_batch_{}_{}", std::process::id(), name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Fixture { dir }
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.join(name)
    }

    fn str_path(&self, name: &str) -> String {
        self.path(name).to_str().unwrap().to_string()
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn batch_matches_per_sample_multi_profile() {
    let fx = Fixture::new("e2e");

    // Two species, two genomes each.
    let (a0, a1) = synth_species(0x9e3779b97f4a7c15);
    let (b0, b1) = synth_species(0x123456789abcdef);
    let db_a = build_species_db(&[a0.clone(), a1.clone()], "a");
    let db_b = build_species_db(&[b0.clone(), b1.clone()], "b");
    // The panels must be big enough for the lowered gates to be meaningful, and each strain
    // must carry private markers, or the comparison below is vacuous.
    for db in [&db_a, &db_b] {
        assert!(db.unique_marker_count(0) >= 10, "strain 0 needs private markers");
        assert!(db.unique_marker_count(1) >= 10, "strain 1 needs private markers");
    }
    let dbs = fx.path("dbs");
    std::fs::create_dir_all(&dbs).unwrap();
    db_a.save(&dbs.join("speciesA.tsv")).unwrap();
    db_b.save(&dbs.join("speciesB.tsv")).unwrap();

    // s1 (paired-end): species A at 70/30 across the mates, species B in R2 only.
    let s1_r1 = fx.path("s1_R1.fq");
    let s1_r2 = fx.path("s1_R2.fq");
    write_fastq(&s1_r1, &synth_reads(&[&a0, &a1], &[4200, 1800], 0xaaa, 150));
    write_fastq(&s1_r2, &synth_reads(&[&b0], &[3000], 0xbbb, 150));
    // The reference for the paired sample: what `cat R1 R2` digests to.
    let s1_cat = fx.path("s1_cat.fq");
    std::fs::write(
        &s1_cat,
        [std::fs::read(&s1_r1).unwrap(), std::fs::read(&s1_r2).unwrap()].concat(),
    )
    .unwrap();

    // s2 (single-end, absolute path in the manifest, no reads2 column).
    let s2_fq = fx.path("s2.fq");
    write_fastq(&s2_fq, &synth_reads(&[&a1, &b0, &b1], &[2000, 1500, 3500], 0xccc, 150));

    let manifest = fx.path("manifest.csv");
    std::fs::write(
        &manifest,
        format!(
            "sample,reads1,reads2\ns1,s1_R1.fq,s1_R2.fq\ns2,{}\n",
            s2_fq.display()
        ),
    )
    .unwrap();

    // batch over both samples at once.
    let batch_out = fx.str_path("batch.tsv");
    let dbs_s = dbs.to_str().unwrap().to_string();
    let manifest_s = manifest.to_str().unwrap().to_string();
    let out = run_strings(&with_gates(vec![
        "batch".into(),
        "--dbs".into(),
        dbs_s.clone(),
        "--manifest".into(),
        manifest_s.clone(),
        "--out".into(),
        batch_out.clone(),
    ]));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("[batch 1/2] sample s1"), "per-sample progress: {stderr}");
    assert!(stderr.contains("[batch 2/2] sample s2"), "per-sample progress: {stderr}");

    let batch_text = String::from_utf8(std::fs::read(fx.path("batch.tsv")).unwrap()).unwrap();
    assert_eq!(
        batch_text.lines().next().unwrap(),
        "#sample\tspecies\tcluster\tabundance\tcoverage\tsupport\tdepth\tglobal_abundance\tsample_fraction\tn_markers"
    );

    // Per-sample reference runs (s1 against the concatenated reads).
    let s1_pred = fx.str_path("s1.pred");
    let s2_pred = fx.str_path("s2.pred");
    let s1_cat_s = s1_cat.to_str().unwrap().to_string();
    let s2_fq_s = s2_fq.to_str().unwrap().to_string();
    for (reads, pred) in [(&s1_cat_s, &s1_pred), (&s2_fq_s, &s2_pred)] {
        run_strings(&with_gates(vec![
            "multi-profile".into(),
            "--dbs".into(),
            dbs_s.clone(),
            "--reads".into(),
            reads.clone(),
            "--out".into(),
            pred.clone(),
        ]));
    }

    // Every batch row, minus the leading sample column, must be byte-identical to the
    // corresponding multi-profile --out row.
    for (sample, pred) in [("s1", &s1_pred), ("s2", &s2_pred)] {
        let from_batch: Vec<String> = batch_text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .filter(|l| l.split('\t').next() == Some(sample))
            .map(|l| l.split_once('\t').unwrap().1.to_string())
            .collect();
        let from_profile = body_lines(Path::new(pred));
        assert!(
            from_batch.len() >= 3,
            "{sample}: expected calls in both species, got {from_batch:?}"
        );
        assert_eq!(
            from_batch, from_profile,
            "{sample}: batch rows must equal multi-profile --out rows field-for-field"
        );
    }
}

/// A manifest row naming a file the reader cannot parse must fail the same way a missing one
/// does. The reader picks FASTA or FASTQ by extension and falls back to FASTA, so a FASTQ
/// called `reads.txt` would otherwise be accumulated as one enormous contig — no error,
/// unbounded memory, and a plausible-looking table at the end.
#[test]
fn batch_unparseable_reads_extension_is_a_hard_error() {
    let fx = Fixture::new("badext");
    let (a0, a1) = synth_species(0x9e3779b97f4a7c15);
    let db_a = build_species_db(&[a0, a1], "a");
    let dbs = fx.path("dbs");
    std::fs::create_dir_all(&dbs).unwrap();
    db_a.save(&dbs.join("speciesA.tsv")).unwrap();

    // The file exists — only its name is unrecognizable.
    let reads = fx.path("s1.txt");
    std::fs::write(&reads, "@r1\nACGT\n+\nIIII\n").unwrap();
    let manifest = fx.path("manifest.csv");
    std::fs::write(&manifest, format!("sample,reads1\nmisnamed,{}\n", reads.display())).unwrap();

    let out = run(&[
        "batch",
        "--dbs",
        dbs.to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--out",
        &fx.str_path("out.tsv"),
        "--enzyme",
        "BcgI",
    ]);
    assert!(!out.status.success(), "an unparseable extension must fail the run");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("misnamed"), "error must name the sample: {stderr}");
    assert!(
        stderr.contains("unrecognized sequence format"),
        "error must say what is wrong: {stderr}"
    );
}

/// A manifest row whose reads file does not exist must abort the whole run with an error
/// naming the sample — never skip the sample and emit a table it is silently absent from.
#[test]
fn batch_missing_reads_is_a_hard_error() {
    let fx = Fixture::new("missing");
    let (a0, a1) = synth_species(0x9e3779b97f4a7c15);
    let db_a = build_species_db(&[a0, a1], "a");
    let dbs = fx.path("dbs");
    std::fs::create_dir_all(&dbs).unwrap();
    db_a.save(&dbs.join("speciesA.tsv")).unwrap();

    let manifest = fx.path("manifest.csv");
    std::fs::write(&manifest, "sample,reads1,reads2\nghost,nope_R1.fq,nope_R2.fq\n").unwrap();

    let out = run(&[
        "batch",
        "--dbs",
        dbs.to_str().unwrap(),
        "--manifest",
        manifest.to_str().unwrap(),
        "--out",
        &fx.str_path("out.tsv"),
        "--enzyme",
        "BcgI",
    ]);
    assert!(!out.status.success(), "missing reads must fail the run");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("ghost"), "error must name the sample: {stderr}");
    assert!(stderr.contains("not found"), "error must say what is wrong: {stderr}");
}
