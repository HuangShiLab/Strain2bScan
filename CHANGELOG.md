# Changelog

All notable changes to Strain2bScan are documented here.

## [Unreleased] — recalibration after the accuracy rewrite

### Fixed — recall regression introduced by the depth fix
- **`--min-abundance` now defaults to 0, not 0.02.** The 0.02 floor was inherited from
  StrainScan and was only survivable while `depth` came from the median over *detected*
  markers, which pins a rare cluster near 1 read/tag no matter how rare it truly is. Removing
  that bias made rare clusters report their real, correctly tiny share — and the unchanged
  floor then deleted them. On a 30x/0.5x mixture the biased estimator reported the rare
  cluster at 3.23% (true 1.64%) and it survived; the unbiased one reports 1.57% and 0.02
  discarded the call outright. This collapsed recall on precisely the samples full of rare
  strains: staggered mocks (MSA1003) and high-host dilutions (99.9% host).

  Precision does not depend on this floor — it is carried by `--min-support`, `--min-coverage`
  and the Layer-1 species gate. Moving it between 0.02 and 0 changes no false positive on the
  mocks, only true positives, and downstream evaluation applies its own presence threshold.
  Pass `--min-abundance 0.02` to restore the previous behaviour.

  The general lesson, now covered by a regression test: **a threshold calibrated against a
  biased estimator has to be recalibrated when the bias is removed.** Fixing the estimator
  alone made the tool measurably more accurate and measurably less sensitive at the same time.

### Fixed — output column order is append-only again
- `profile` writes `cluster/abundance/coverage/support/depth/n_markers` and `multi-profile`
  writes `species/cluster/abundance/coverage/support/depth/global_abundance/sample_fraction/
  n_markers`. The previous rewrite inserted `depth` as the third column, shifting `coverage`
  and `support` by one; positional readers do not fail on that, they just silently produce
  wrong numbers. The leading columns now match the pre-rewrite layout exactly and every new
  quantity is appended.

## [Unreleased] — accuracy and performance rewrite

Accuracy first, then speed. Existing databases remain readable — marker values are still
FNV-1a of the canonical tag, and the on-disk format is unchanged.

### Added — `sample_fraction`, a cross-sample-comparable abundance
- `multi-profile` gains a `sample_fraction` column: the cluster's share of **all** tag
  observations in the sample, `depth × n_markers / Σ counts`. `global_abundance`'s denominator
  is whatever that run resolved, so it is not comparable between samples — change the depth,
  the reference panel, or how much of the community has no reference, and the denominator
  moves. In a defined mock that is harmless; in a real sample, where most of the community may
  be unreferenced, it silently rescales everything. `sample_fraction`'s denominator is fixed by
  the sequencing instead, and the remainder is reported rather than hidden:

      coverage of sample: strain calls account for 34.2% of 1234567 tag observations
      (65.8% unclassified — unresolved species, no reference, host, error)

- The two columns answer **different biological questions**, and both are now available. In the
  standard metagenomics terminology (CAMI, and the ATCC microbiome standards):

  | ground truth is stated as | meaning | compare against |
  |---|---|---|
  | **taxonomic abundance**, cell/organism abundance | fraction of *cells* | **`global_abundance`** |
  | sequence abundance, genomic/DNA abundance | fraction of *DNA* | `sample_fraction` |

  `global_abundance ∝ depth` is a cell fraction because unique markers are **single-copy**: a
  genome with more tags yields proportionally more reads, so reads-per-tag cancels the genome
  size out (`depth = reads/G ∝ cells·G/G`). `sample_fraction ∝ depth × n_markers` puts the
  genome size back in and is therefore a DNA fraction. On a mock with two clusters at equal cell
  depth but 2.0 Mb and 1.6 Mb genomes, the columns correctly report 0.510/0.490 (cells) and
  0.564/0.436 (DNA). ATCC MSA-1002 states **taxonomic abundance**, so benchmark it against
  `global_abundance`; using the wrong column imposes a systematic error scaling with the
  genome-size spread of the panel (>6× across MSA-1002).

### Fixed — quantification
- **`multi-profile` now uses the cross-species evidence for quantification, not just gating.**
  Cluster-uniqueness is only defined *within* one species database, so a tag can be unique to a
  cluster there and still occur in a congener's genomes — and when that congener is co-present,
  its reads land on the tag and inflate the cluster's depth. Panels routinely contain such pairs
  (*S. aureus*/*S. epidermidis*, the three streptococci and the two lactobacilli in ATCC
  MSA-1002, and most oral communities), so this was a systematic abundance error rather than a
  rare accident. Detection and depth are now restricted to markers that are specific to their
  species across the whole panel; the restriction is automatic and reported in the run header.
  On a mock with two co-present congeners sharing a region, the affected cluster's depth was
  overstated **3×** (29.9× vs. a true 10×); with the filter it is 10.2×. Bray–Curtis similarity
  to truth: global 0.791 → **0.996**, within-species 0.936 → **0.997**.
  `--no-cross-species-filter` restores the old behaviour for comparison.
- **Absolute depth is now recoverable, so a cross-species composition is possible at all.**
  Previously each species DB was profiled independently and the only number reported was a
  within-species fraction, so a species at 10× the depth of another reported identical
  abundances and all between-species information was discarded — a community composition simply
  could not be reconstructed from the output. `multi-profile` now reports the absolute per-tag
  `depth` plus a `global_abundance` column derived from it, alongside the within-species
  `abundance` (which remains the primary column). On a 2-species/4-cluster mock with a 10×
  species-level split, Bray–Curtis similarity of the global composition to truth is **0.991**;
  reconstructing it the only way the old output allowed (concatenate + renormalize) gives
  **0.591** (L1 error 0.017 vs 0.818).
- **Depth estimator no longer flattens rare strains.** `unique_marker_depth` took the median
  over *detected* markers only, so a cluster at 0.05× and one at 1× both reported a depth of
  ~1. It is now the **zero-inclusive** mean over the whole unique-marker panel, trimmed at the
  top 1% to retain the median's robustness to collapsed repeats.

### Fixed — sensitivity
- **The singleton filter and the species floor scale with depth.** `count >= 2` and the
  200-marker species floor were fixed constants, both unreachable at low depth: at 0.05× only
  ~5% of any marker panel is observable at all, so a genuinely present species was filed as
  unresolvable no matter how clean the data was. Both now scale with the reachable fraction
  `1 − e^(−λ)`. On a 5% subsample of the mock (~0.1–2× depth) recall went from **2/4 to 3/4
  clusters**, with the global composition at Bray–Curtis similarity 0.970 to truth.
  `--fixed-gate` restores the previous behaviour.

  Two bounds keep this from trading away precision:
  - The species floor relaxes only to `MIN_FLOOR_FRACTION` (25%) of its configured value. An
    unbounded scaling is self-cancelling — the observed count is itself proportional to the
    reachable fraction, so `present >= floor·r` reduces to `total >= floor` at every depth,
    leaving `--min-species-detect` as the only real threshold.
  - The **coverage floor is deliberately NOT depth-scaled.** Scaling it is provably inert:
    depth is estimated from the same panel whose breadth is being tested, and whenever every
    count is 0 or 1 the scaled gate is satisfied for any coverage including zero. It stays
    absolute and remains the guard that stops a large, sparsely-hit panel from being called;
    lower `--min-coverage` to trade precision for recall.

### Added
- **Streaming FASTA/FASTQ with gzip support.** `.fq.gz` / `.fna.gz` are read directly
  (piped through `gzip -dc`, so the crate stays dependency-free). Samples are digested in
  batches instead of being loaded whole: on a 289 MB / 4M-read sample peak RSS fell from
  **282 MB to 12 MB**, and gzipped input costs no extra wall-clock.
- `--min-abundance` and `--fixed-gate` flags; `--out` for `multi-profile`.
- `depth` column in all prediction output (absolute reads-per-tag, comparable across species).

### Changed
- **Output format** — `profile` TSV is now `cluster/abundance/depth/coverage/support`, and
  `multi-profile` writes `species/cluster/abundance/global_abundance/depth/coverage/support`,
  grouped by species with the most abundant cluster first. `abundance` keeps its original
  within-species meaning; `global_abundance` and `depth` are new. Downstream plotting scripts
  need updating (and any that concatenated per-species fractions should read
  `global_abundance` instead of re-deriving).
- `--min-abundance` applies **within** a species, matching the primary column. Whether a
  species is present at all remains the Layer-1 species gate's decision.
- Enzyme matching is case-insensitive, so soft-masked reference genomes digest correctly.
- `Params` lost the dead `alpha` / `l1_ratio` / `max_iter` fields (never read; the Elastic Net
  solver takes them as arguments and remains available but off the main path).
- Internals: FxHash replaces SipHash for `u64`-keyed maps; tag canonicalization and digestion
  are allocation-free. Combined with streaming, profiling 4M reads is **4.7× faster** with one
  enzyme and **3.4× faster** with all 16.

## [0.1.0] — 2026-06-30

Initial release: a Rust strain-level metagenomic profiler that applies the StrainScan
resolution framework to 2bRAD tag markers.

### Added
- **Type-IIB digestion** for all 16 enzymes (BcgI, CspCI, AlfI, …); single-enzyme (BcgI
  2bRAD data) and multi-enzyme (`--enzyme all`, digital digestion of shotgun reads) modes.
- **Single-copy 2b-tag markers**; canonical `u64` encoding; sparse strain×marker database
  with a unique-marker inverted index.
- **Within-species clustering (CST):** single-linkage at 0.95 similarity; exact Jaccard for
  small panels and bottom-k **MinHash** sketches for large panels (identical partitions on
  real data, near-linear build).
- **Marker classification:** species-core / shared-partial / cluster-specific / strain-specific.
- **Layer-2 profiling:** unique-marker presence detection + non-negative Elastic Net
  abundance estimation; depth from a strain's unique markers.
- **Multi-species mode (`multi-profile`):** digest a sample once, match all per-species DBs
  in parallel; a Layer-1 species gate (species-specific markers) suppresses cross-species
  false positives.
- **Parallelism:** dependency-free multi-threaded digestion/clustering (`STRAIN2BSCAN_THREADS`).
- **Evaluation:** precision/recall/F1, L1 and Bray–Curtis via the `evaluate` subcommand.
- CLI: `build`, `cluster`, `profile`, `multi-profile`, `info`, `evaluate`, `demo`, `cst-demo`.
