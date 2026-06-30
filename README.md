# Strain2bScan

**Fast strain-level metagenomic profiling on 2bRAD-reduced k-mer markers** — a Rust
reimplementation of the [StrainScan](https://github.com/liaoherui/StrainScan) resolution
framework (clustering + unique-marker Layer-2) operating on **2bRAD tags** instead of the
full k-mer set. It trades a controlled amount of low-depth sensitivity for large gains in
speed and memory, and adds native handling of BcgI 2bRAD experimental data.

> Built on the ideas of StrainScan (Liao et al., *Microbiome* 2023) and the 2bRAD tag
> extraction of [Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M).

## Why 2bRAD markers

StrainScan resolves strains by scoring **unique markers** (k-mers specific to a strain or
cluster). 2bRAD type-IIB digestion yields a sparse, reproducible ~1–2% subset of the genome —
exactly the kind of low-redundancy marker set that algorithm wants. Using tags instead of all
k-mers shrinks the database ~50–100×, so digestion and matching are far faster and lighter,
while clustering and Layer-2 logic are preserved.

## Highlights

- **All 16 type-IIB enzymes**; single-enzyme (BcgI 2bRAD data) or multi-enzyme digital
  digestion of shotgun reads (`--enzyme all`).
- **Within-species clustering (CST):** single-linkage at 0.95; exact Jaccard for small panels,
  **MinHash sketches** for large ones (identical partitions on real data, near-linear build).
- **Layer-2 profiling:** present-cluster detection from unique markers + robust abundance from
  unique-marker depth (non-negative Elastic Net solver also included).
- **Multi-species mode (`multi-profile`):** digest a sample **once**, match all per-species
  DBs in parallel; a Layer-1 **species gate** suppresses cross-species false positives.
- **Parallel** (dependency-free, `std` threads) — ~10× build, ~7× profile on 16 cores.
- Honest reporting: tells you when a species is detectable but **not strain-resolvable** with
  the given enzyme set (e.g. BcgI alone on a low-diversity species).

## Install

```bash
git clone https://github.com/HuangShiLab/Strain2bScan
cd Strain2bScan
cargo build --release      # binary at target/release/strain2bscan
cargo test                 # 23 tests
```

## Two input modes (`--enzyme`)

`all` (16 enzymes), one enzyme (`BcgI`), or a list (`BcgI,CspCI`). The genome DB and the
sample must use the same set — `profile` reads it from the DB header automatically.

1. **BcgI 2bRAD data** (reads already are 2b tags): `--enzyme BcgI`. Sparse; if a species
   lacks cluster-specific tags, `cluster` reports `✗ NOT DOABLE` and `profile` reports the
   species as detectable but not strain-resolvable.
2. **Conventional metagenome** (150 bp / long reads): `--enzyme all` — digitally digest with
   all 16 enzymes to enrich strain markers (~hundreds× more tags).

## Usage

```bash
# single species: build a cluster DB, profile a sample, evaluate
strain2bscan cluster      --genomes acnes_genomes/ --enzyme all --out acnes.db.tsv --similarity 0.95
strain2bscan profile      --db acnes.db.tsv --reads sample.fq --out pred.tsv
strain2bscan evaluate     --pred pred.tsv --truth truth.clusters.tsv --present 0.01

# many species at once: one DB per species in a dir; sample digested once, matched in parallel
strain2bscan multi-profile --dbs species_dbs/ --reads sample.fq --enzyme all --min-species-markers 200

# self-contained demos (no data needed)
strain2bscan demo          # conspecific 70/30 mixture resolved by Layer-2
strain2bscan cst-demo      # 4 genomes → 2 clusters; marker classes; cluster profiling

STRAIN2BSCAN_THREADS=8 strain2bscan cluster ...   # control threads (default: all cores)
```

`cluster` also writes `<out>.members.tsv` (genome→cluster) for remapping ground truth.

## Marker taxonomy: strain-specific vs species-specific

Within a species each tag is classified by its **within-species incidence**: `SpeciesCore`
(in all clusters → detects the species, not strains), `ClusterSpecific` (one cluster),
`StrainSpecific` (one genome), `SharedPartial`. Cluster/strain-specific tags are the Layer-2
markers. They are derived from **all** 2b tags of the species' genomes (StrainScan's
all-k-mer approach) — *not* from a species-unique marker set. Fast2bRAD-M's species-unique
markers are computed by comparing each genome against genomes of *other* species (for species
detection) and are orthogonal to within-species strain structure.

## Modules

| file | role |
|---|---|
| `enzymes.rs` | type-IIB digestion (all 16 enzymes) + enzyme-set parsing |
| `markers.rs` | canonical tag → `u64` marker; single-copy filter; (parallel) FASTA/FASTQ digest |
| `db.rs` | sparse strain×marker DB with unique-marker index |
| `cst.rs` | within-species clustering (exact + MinHash) and marker classification |
| `identify.rs` | Layer-2: unique-marker detection + unique-marker-depth abundance (+ NNLS) |
| `parallel.rs` | dependency-free parallel map (`std::thread::scope`) |
| `bench.rs` | precision/recall/F1, L1, Bray–Curtis metrics |
| `main.rs` | CLI |

## Performance (vs StrainScan, real C. acnes)

| | Strain2bScan | StrainScan |
|---|---|---|
| profile / sample | **~3.6 s (0.5 s @16 threads)** | ~6.8 s |
| peak memory | **~105 MB** | ~830 MB |

Multi-species cost is independent of species count (digest once, match many); per-sample
throughput is linear. Full benchmarks, scaling curves and accuracy: see the
**[Strain2bScan-paper](https://github.com/HuangShiLab/Strain2bScan-paper)** repository.

## Scope & honest caveats

- 2bRAD tags are ~50–100× sparser than full k-mers, so Strain2bScan is **less sensitive at
  <1× per-strain depth** than full-k-mer StrainScan; it matches accuracy at sufficient depth
  (≈≥5×) while being much faster/lighter. The wet-lab-2bRAD low-input advantage is a separate
  use case.
- Reads are matched exactly; an error-tolerant mode and gzip/streaming I/O
  ([needletail](https://crates.io/crates/needletail)) are planned.
- Species selection (which species to resolve at strain level) comes from Fast2bRAD-M's
  species-level profiling output; Strain2bScan then digests those species' genomes itself.

## Citation

If you use Strain2bScan, please cite this repository and StrainScan
(Liao et al., High-resolution strain-level microbiome composition analysis from short reads,
*Microbiome* 2023; doi:10.1186/s40168-023-01615-w).

## License

MIT — see [LICENSE](LICENSE). Strain2bScan reimplements StrainScan's framework (MIT) on
2bRAD tags; StrainScan and Fast2bRAD-M are credited above.
