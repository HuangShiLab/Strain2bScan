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
- **Layer-2 profiling:** present-cluster detection from unique markers + **absolute** per-tag
  depth from a zero-inclusive trimmed mean over the unique-marker panel (non-negative Elastic
  Net solver also included, but not on the main path).
- **Multi-species mode (`multi-profile`):** digest a sample **once**, match all per-species DBs
  in parallel. A Layer-1 **species gate** suppresses absent species, and detection/quantification
  are automatically restricted to markers specific to their species *across the whole panel* —
  without that, a co-present congener's reads land on tags that merely look cluster-specific
  within one database and inflate its depth (3× on a two-congener mock). Abundance is
  reported at **three** scopes, because per-species fractions cannot simply be concatenated into
  a community composition:
  | column | denominator | question it answers |
  |---|---|---|
  | `abundance` | this species | given the species is present, how does it split across strains? (**primary**) |
  | `global_abundance` | what this run resolved | composition of the resolved part — fine for a defined mock, **not comparable between samples** |
  | `sample_fraction` | all tag observations | share of the sequencing; denominator fixed by the data, so it **is** comparable between samples |

  `global_abundance ∝ depth` is a **cell** fraction (what CAMI and the ATCC standards call
  *taxonomic abundance*) — unique markers are single-copy, so reads-per-tag cancels genome size
  out. `sample_fraction ∝ depth × n_markers` puts it back in and is a **DNA** fraction (*sequence
  / genomic abundance*). Match the column to how your ground truth is stated: ATCC MSA mixes give
  taxonomic abundance, so benchmark against `global_abundance`. The unclassified remainder is
  printed rather than hidden.
- **Depth-adaptive gating:** the singleton filter and the species-marker floor scale with the
  estimated per-tag depth (`1 − e^(−λ)`), so low-input and high-host samples are not gated out
  by thresholds unreachable at their depth. The relaxation is bounded, and the coverage floor
  stays absolute, so precision is not traded away. `--fixed-gate` restores the old constants.
- **Streaming I/O with gzip:** `.fq.gz` / `.fna.gz` read directly; peak memory is one batch,
  not the whole file (still dependency-free — decompression pipes through `gzip -dc`).
- **Assembly-quality aware:** variable completeness biases Jaccard toward *spurious splits*
  (an incomplete genome looks distant from its complete twin). `build`/`cluster` always flag
  likely-incomplete genomes (single-copy tag count ≪ the conspecific median) and can drop
  fragmented/incomplete assemblies before clustering (`--max-contigs N`, `--min-tag-fraction F`).
- **Parallel** (dependency-free, `std` threads) — ~10× build, ~7× profile on 16 cores.
- Honest reporting: tells you when a species is detectable but **not strain-resolvable** with
  the given enzyme set (e.g. BcgI alone on a low-diversity species).

## Install

```bash
git clone https://github.com/HuangShiLab/Strain2bScan
cd Strain2bScan
cargo build --release      # binary at target/release/strain2bscan
cargo test                 # 47 tests
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

## Abundance: exact definitions

Four quantities are reported and they answer different questions. Mixing them up is the
easiest way to get a wrong benchmark, so they are defined here exactly as implemented.

### Notation

| symbol | meaning |
|---|---|
| `j` | a cluster (the strain unit); `s(j)` is its species |
| `M_j` | **all** single-copy 2b tags carried by cluster `j`; `G_j = \|M_j\|` (reported as `n_markers`) |
| `U_j` | cluster `j`'s **discriminating** markers: tags in `M_j` carried by no other cluster and — under `multi-profile` — by no other species. `N_j = \|U_j\|` |
| `c_m` | how many times tag `m` was observed in the sample |
| `T` | `Σ_m c_m` over every tag observed in the sample (the whole digest, classified or not) |

Only single-copy tags enter `M_j` (`markers::single_copy_markers`), so one genome copy
contributes exactly one copy of each tag. That is what makes depth proportional to genome
copy number.

### Per-cluster primitives

    D_j  = |{ m ∈ U_j : c_m ≥ 1 }|                       detected markers
    κ_j  = the (⌊D_j/100⌋ + 1)-th largest of { c_m }_{m ∈ U_j}     winsorization cap

    depth_j    = (1 / N_j) · Σ_{m ∈ U_j} min(c_m, κ_j)            reads per tag
    coverage_j = D_j / N_j                                        breadth

`depth_j` is a mean over the **whole** panel `U_j`, zeros included — a marker that was not
observed is evidence the cluster is rare, not a missing data point. Averaging only over
detected markers (or taking their median) makes a cluster at 0.05x and one at 1x both report
≈1, which flattens the entire composition. The cap `κ_j` winsorizes the top 1% of *detected*
observations so a collapsed repeat cannot inflate the mean; it caps rather than discards,
because discarding `N_j/100` entries deletes real signal once breadth falls below ~1%.

Note `depth_j` is estimated on `U_j` (discriminating markers only, so the reads are
unambiguously this cluster's), but the cluster physically carries all of `M_j`. Both sets are
used below, and using the wrong one is the subtle error to avoid.

### The four reported columns

**1. `abundance` — within-species (primary)**

    abundance_j = depth_j / Σ_{k ∈ s(j)} depth_k

over the clusters of species `s(j)` that passed the support and coverage gates, then
low-abundance clusters are dropped (`--min-abundance`, **default 0**) and the survivors
renormalized. Sums
to 1.0 within each species. Answers: *given this species is present, how does it split across
strains?* The denominator is the species itself, so it is comparable between samples.

**2. `species_abundance` / `global_abundance` — community level**

    species_abundance_s = (Σ_{k ∈ s} depth_k) / Σ_k depth_k
    global_abundance_j  = depth_j / Σ_k depth_k        (k over every emitted cluster)

and these compose exactly:

    global_abundance_j = species_abundance_{s(j)} × abundance_j

so the flat computation and the hierarchical one (species abundance first, then split within
species) give the same number — verified to 1e-6 on the mock. Substituting an externally
computed species abundance (e.g. Fast2bRAD-M's) into the identity is therefore valid.

⚠️ The denominator spans **only the clusters actually emitted**. Species that were detected
but not strain-resolvable, clusters below the gates, organisms with no reference, host DNA
and sequencing error are all outside it. In a defined mock (membership known, everything
detectable) that is harmless; in a real sample the denominator changes with depth and with
how much of the community is unreferenced, so `global_abundance` is **not comparable between
samples**.

**3. `sample_fraction` — share of the sequencing**

    sample_fraction_j = depth_j · G_j / T

The denominator is fixed by the sequencing rather than by how well profiling went, so this
one **is** comparable between samples, and the shortfall from 1.0 is reported explicitly as
unclassified. The numerator converts a per-tag depth into a tag mass: `depth_j` alone is
reads-per-tag and cannot be compared against a read total. Clusters sharing a tag each
contribute their own depth to it, so the shares add without double counting — on a closed
mock where every organism has a reference, the column sums to 100.0%.

### Cells vs DNA — pick the column your ground truth uses

`depth ∝ genome copies`, and `depth · G ∝ genome copies × genome size ∝ DNA mass`. So:

| column | proportional to | biological meaning |
|---|---|---|
| `abundance`, `global_abundance` | `depth` | **taxonomic / cell abundance** (genome copies) |
| `sample_fraction` | `depth · G` | **DNA / sequence mass** |

They differ by genome size. On a mock with two clusters at equal cell depth but 2.0 Mb and
1.6 Mb genomes, the columns correctly report 0.510/0.490 and 0.564/0.436. A ground truth
stated as *taxonomic abundance* (ATCC MSA mixes) must be compared against `global_abundance`;
one stated as equal *genomic DNA* against `sample_fraction`. Using the wrong one imposes a
systematic error that scales with the genome-size spread of the panel — for a 20-species mock
spanning 0.7–4.6 Mb that is a factor of ~6 between the extremes.

### Assumptions

- Uniform tag recovery across loci. Real 2bRAD has enzyme/GC bias, which inflates the
  variance of `depth_j` but not (to first order) its expectation.
- `G_j` counts single-copy tags only, so `depth · G` under-estimates true tag mass by the
  multi-copy fraction. This is systematic and roughly constant across organisms, so it
  largely cancels in `sample_fraction` ratios.
- `depth_j` is measured on `U_j` and extrapolated to `M_j`. Valid when discriminating markers
  are not systematically biased in coverage relative to the rest of the genome.

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
| `quality.rs` | assembly-quality filter (contig count + tag-count completeness proxy) |
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
- Reads are matched **exactly** — an error-tolerant digestion mode is still planned. (Streaming
  and gzip I/O are implemented; no `needletail` dependency was needed.)
- Species selection (which species to resolve at strain level) comes from Fast2bRAD-M's
  species-level profiling output; Strain2bScan then digests those species' genomes itself.

## Citation

If you use Strain2bScan, please cite this repository and StrainScan
(Liao et al., High-resolution strain-level microbiome composition analysis from short reads,
*Microbiome* 2023; doi:10.1186/s40168-023-01615-w).

## License

MIT — see [LICENSE](LICENSE). Strain2bScan reimplements StrainScan's framework (MIT) on
2bRAD tags; StrainScan and Fast2bRAD-M are credited above.
