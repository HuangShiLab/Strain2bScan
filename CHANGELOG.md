# Changelog

All notable changes to Strain2bScan are documented here.

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
