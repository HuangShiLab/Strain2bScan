# Changelog

All notable changes to Strain2bScan are documented here.

## [Unreleased] — the tree descent was bypassing the cross-species filter

### Fixed — `--layer1 cst` scored unfiltered markers under `multi-profile`

`multi-profile` restricts detection and quantification to the markers specific to each
species across the whole panel, so a co-present congener's reads cannot land on a tag that
merely looks cluster-specific inside one database. `restrict_to` records that as a
`quant_mask` on the database and `unique_markers` consults it — but the tree's node sets are
reached through `Cst`, which the mask never touched. `descend_tree` does not even receive the
database. So the flat path was filtered and the tree path was not, and the two were not
scoring the same sample: every internal node was tested on markers this species shares with
its congeners, which is exactly the cross-talk the filter exists to stop.

`descend_tree_masked` takes the mask and applies it to every node set it reads — the fire
test, the informativeness test, the pooled panel, and the clade fallback. `tree_utility`
counts informative nodes the same way, so `auto` decides on the panel the descent will
actually score rather than on a larger one.

Single-species `profile` carries no mask and is byte-identical on all five C. acnes samples
across `auto`, `unique` and `cst`. Only `multi-profile` changes — and any measurement of
`--layer1 cst` under `multi-profile` taken before this needs re-running, because the tree had
access to evidence the flat path was denied.

### Added — `multi-profile` reports the Layer-1 decision per species

`auto` decides per species, so one global line could not report it; the only way to tell
which species took the tree was to diff whole runs. Each resolved species now prints its
decision and the two counts behind it.

## [Unreleased] — make Layer-2 reachable, choose Layer-1 from the database

### Fixed — `--layer2 enet` could never reach the case it was written for

The joint fit exists to resolve a cluster whose tag set is contained in a co-present
relative's: it has no unique markers, so the flat path cannot see it, but the shared rows
observe `w_A + w_B` while the A-only rows observe `w_A`, and the two together give `w_B`.
That capability was unreachable from the CLI. `profile()` handed the fit only the clusters
Layer-1 had already called, and `profile_unique` skips any cluster with an empty unique panel
(`panel == 0`) — precisely the clusters in question. With one column, `calls.len() > 1` was
false and the layer did not even run.

The layer's own test passed because it hands `build_l2_design` the candidate list `[0, 1]`
directly, stepping over the wiring it was meant to exercise. A second test now runs the same
scenario through `profile()`, which is what the CLI calls.

Two changes make it reachable. `subset_candidates` collects the clusters Layer-1 structurally
cannot see and offers them to the fit alongside the called ones. And those columns bypass the
pre-scan: it ranks by *residual markers explained* and consumes the winner's markers, so a
cluster contained in one already chosen scores zero on the next iteration and can never be
selected at any threshold — breadth is the wrong instrument for a candidate distinguished by
depth. The non-negative solve is free to give them zero weight instead.

Admission is deliberately narrow, since a cluster with no unique evidence is one whose every
supporting read is also a relative's:

- it must have too few unique markers to have been judged on its own, so this can never
  overturn a Layer-1 rejection made on unambiguous evidence;
- it must sit inside a **single** called cluster to within fewer markers than the support
  floor. A *fractional* containment cannot express this on a within-species panel — two
  conspecific genomes already share well over 95 % of their tags, so a 0.95 threshold is met
  by ordinary relatedness. Measured here, the fractional form admitted the same spurious
  cluster in 4 of 5 samples (precision 1.000 → 0.793); the absolute form admits none;
- the fit must give it at least `MIN_SUBSET_SHARE` of the total depth.

### Changed — `--layer1` defaults to `auto`, decided per database

Whether a Cluster Search Tree helps is a property of the panel, not of the software, so a
global flag is the wrong shape for it. `auto` reads the answer off the database: descend only
if some cluster falls below the support floor — the only case pooling can change an outcome —
**and** some internal node carries enough markers to pool. Both halves are necessary, and on a
dense panel the second fails: of 542 internal nodes on 543 *C. acnes* genomes, 373 carry zero
group-specific markers, because clustering at τ already merges anything similar enough for a
clade to have a distinct core.

The decision and the two counts behind it are printed on every run; an automatically chosen
algorithm must not be a silent one. `--layer1 unique|cst` still forces either path. A database
with no tree resolves to `unique`, so nothing changes for databases built by `build`.

### Note — with Layer-2 reachable, it is measurably harmful here

Previously `--layer2 enet` looked inert because detection metrics cannot see it — it never
changes *which* clusters are called. Scored on abundance accuracy over the clusters both
agree on, it is not inert but wrong: Bray–Curtis 0.035 → 0.127 and mean absolute relative
error 0.183 → 0.794 against the default estimator, on identical detections.

The cause is collinearity, and it is structural rather than a tuning problem. On this panel
each cluster carries ~33 100 markers of which only 29–115 are unique, so the design columns
are ~99.7 % identical; the shared rows constrain the *sum* of two near-identical clusters and
say almost nothing about the split. Penalising does not help — sweeping `--enet-alpha` makes
it monotonically worse (BC 0.127 at 0, 0.145 at 0.001, 0.223 at 0.01, 0.330 at 0.1, 0.360 at
1.0), which agrees with the earlier finding that 0 is the best setting. The few unique markers
carry the split directly and the many shared ones do not, so the unique-marker mean wins
despite using 0.3 % of the data.

`depth` therefore remains the default. The layer is now reachable and honestly measured rather
than unreachable and assumed fine.

## [Unreleased] — lower the support floor to 8

### Changed — `min_support_markers` 10 → 8

This gate sets the detection limit, and 10 was a round number rather than a measured one.
Support tracks `panel × (1 − e^(−λ))` almost exactly, and on 2bRAD markers both factors are
small: the discriminating panel is a few dozen tags (median 53 across a 419-cluster
*C. acnes* panel, 13–115 for the clusters actually present), and a strain at 5 % of a sample
sequenced to ~5× per tag sits at λ ≈ 0.27, where only ~24 % of any panel is observable. So
~35 markers × 24 % ≈ 8 observed — and a floor of 10 discards strains that are plainly there.

Measured on that panel, three of 15 truth clusters were lost on this gate alone, all of them
at 2.7–6.3 % relative abundance, and none of them rejected by coverage or consistency:

| cluster | truth | panel | depth | coverage | support |
|---|---|---|---|---|---|
| C249 | 2.7 % | 40 | 0.200× | 20.0 % | 8 |
| C229 | 6.3 % | 34 | 0.294× | 23.5 % | 8 |
| C202 | 5.1 % | 28 | 0.143× | 14.3 % | 4 |

Sweeping the floor on the flat path:

| floor | precision | recall | AUPR |
|---|---|---|---|
| 10 | 1.000 | 0.800 | 0.800 |
| **8** | **1.000** | **0.950** | **0.950** |
| 6 | 0.900 | 0.950 | 0.950 |
| 4 | 0.780 | 1.000 | 0.983 |

8 is the last free step: it recovers C249 and C229 with precision unchanged at 1.000, while 6
and 4 start trading precision away. C202, at support 4, needs a floor of 4 and is where this
panel and depth genuinely run out of signal. Bray–Curtis improves 0.058 → 0.042.

Calibrate this on your own panel — the argument above is arithmetic, but the number 8 comes
from one species at one depth. `--min-support N` overrides it per run.

Worth noting what this floor does at the other end: C204, the *dominant* strain of sample 5 at
83.8 % abundance, 92.3 % coverage and 4.2× depth, cleared the old floor with exactly 10
markers because its panel is only 13 — one dropout from being lost — while a cluster with
8 000 markers needs the same 10, which is 0.1 % of its panel. A floor that scales with panel
size would separate those two regimes; that is a larger change and is not made here.

## [Unreleased] — fix three silent failures in the Cluster Search Tree descent

Measured on the first real panel the port has seen: 543 *Cutibacterium acnes* genomes
(`MockMetagenomes4Benchmark/Cutibacterium_acnes_0.01`, 5 samples, 100 k read pairs, 2–5
genomes each), `--enzyme all`. Ground truth is stated per genome and remapped through
`<out>.members.tsv`. Because the 14 truth genomes merge at Jaccard 0.9825–0.9961 they all
land in one cluster at the default 0.95, so accuracy is reported at `--similarity 0.997`,
where they resolve to 13 clusters.

`--layer1 cst` scored **below** the flat path it is supposed to extend — precision 0.800,
recall 0.600, AUPR 0.600 against the flat path's 1.000 / 0.800 / 0.800. Three separate
causes, each fixed here. Afterwards CST matches the flat path exactly on every sample and at
every support floor; clustering is unchanged and all 53 tests still pass.

### Fixed — a node with exactly `MIN_NODE_MARKERS` markers silently pruned a present subtree

`MIN_NODE_MARKERS` was 10, deliberately set equal to `min_support_markers` "so the two gates
stay on one scale". That is the bug. A node is allowed to veto its whole subtree once it is
*informative*, but it only fires when `support >= min_support_markers` — so a node holding
exactly 10 markers was trusted to prune while being unable to fire unless all 10 were seen at
count ≥ 2. One Poisson dropout cut a present subtree, with no escape hatch, because
`informative` was true and the "cannot be tested, descend anyway" branch did not apply.

Observed on sample 1: node N831 held 10 markers, 9 detected, coverage 0.90 at depth 4.7× —
overwhelming evidence of presence — and failed by one marker, cutting the 331-leaf subtree
that held the answer. The descent then reported a clade spanning 332 of 419 clusters.

`MIN_NODE_MARKERS` is now 30, and documented as needing to sit well **above** the support
floor rather than on it: at a floor of 10 a 30-marker node still clears it with a third of its
panel missing. Nodes below 30 are treated as untestable and entered, which costs work but
cannot lose a true subtree.

### Fixed — a leaf's tree panel was the intersection while the flat path scores the union

`cluster_db` gives a cluster the **union** of its members' markers; `build_tree` gave the same
cluster the **intersection**. So `--layer1 cst` scored every multi-genome cluster on a strictly
smaller panel than the flat path uses, and for a cluster whose members disagree enough the
intersection empties out and the cluster becomes uncallable at any threshold — C247 (5
members) had 0 markers by intersection against 115 by union, and was missed on sample 4 while
the flat path found it without difficulty.

Leaves now use the union, so the tree is a strict extension of the flat path rather than a
weaker competitor. Internal nodes keep the intersection: a marker is diagnostic of a clade
only if every genome under it carries it. **Databases must be rebuilt** — the leaf marker sets
are serialized.

### Fixed — the clade fallback climbed to the root

Nothing bounded the upward walk, so when a region of the tree was rejected it climbed to a
near-root node whose panel is the species core — fully covered in any sample of that species,
so it passes every gate — and reported that. New `MAX_FALLBACK_CLADE` (8) stops the climb once
the clade is too broad to be an answer; above it, report nothing and let the species layer say
"present, not resolvable", which is the truth. Inert once the first fix is in, kept as a guard.

### Note — the tree still does not beat the flat path on this panel, and `diagnose-tree` said so

After all three fixes CST equals the flat path and never exceeds it, at support floors 10, 8,
6 and 4 alike. `diagnose-tree` predicted this before the benchmark ran: of 542 internal nodes,
373 carry **zero** group-specific markers, and the verdict is `NOT viable (median 0 < 25)`.
The nodes that do carry markers merge at similarity ≈ 1.0 — near-duplicates that collapse into
one cluster anyway — so the tree has signal only where it is not needed. Relaxing the clade
definition does not rescue it: redefining an internal node's markers as a prevalence core
(present in ≥ 90 % of descendants rather than all) moves the median from 0 to 2, against a
floor of 25. With genomes this similar a clade-shared marker is almost always carried outside
the clade too, and the "minus everything outside" subtraction empties the node.

`--layer2 enet` never changed a detection — it only runs when Layer-1 returns more than one
candidate and only refits depths — and it made abundance worse (Bray–Curtis 0.058 → 0.101;
sample 4 0.100 → 0.296). It is left alone here.

## [Unreleased] — expose the recommended enzyme set

### Added — `--enzyme recommended`
- The 14 non-degenerate enzymes: everything in the table except **HaeIV** and **Hin4I**.

  `enzymes.rs` has always documented these two as "highly degenerate and excluded from the
  recommended multi-enzyme set (noisy)", but that set was never selectable, so multi-enzyme runs
  had to use `all` and take them along. Measured on a 1.8 Mb genome:

  | set | tags | 1 per | cluster-unique |
  |---|---|---|---|
  | `BcgI` | 984 | 1 829 bp | 14.5% |
  | `recommended` (14) | 16 016 | 1 116 bp | 13.9% |
  | `all` (16) | 23 815 | 1 76 bp | 14.0% |

  HaeIV and Hin4I yield ~4 150 tags each — comparable to the most productive specified enzymes —
  because their IUPAC-degenerate sites match far more often. Together they are 35% of `all`'s
  panel. That looseness is also what makes their sites the easiest for a sequencing error to
  create or destroy, which adds noise to the depth estimate and to the between-strain
  discrimination that depends on it. The unique-marker *fraction* is flat (~14%) across all three
  sets, so dropping them costs a third of the panel without changing its composition.

  Whether that trade is net-positive is an empirical question per dataset, and now it is one flag
  to test rather than 14 enzyme names to type. Note that marker count is unlikely to be the
  binding constraint on the current mock precision: a shadow cluster's overlap fraction `f` is a
  property of the genomes, not of how densely they are sampled, so denser digestion does not
  shrink it.

## [Unreleased] — reject shadow clusters (depth–breadth consistency)

### Fixed — the dominant same-species false positive
- **`--min-consistency` (default 0.5)** rejects a cluster whose observed markers are too *deep*
  for how *few* of them were seen: `coverage / (1 − e^(−depth))`.

  Breaking the per-sample false positives down by species showed two unrelated populations, and
  only one of them was trace contamination. The other — the dominant one — is a **shadow**: the
  strain in the sample is not exactly any reference, so it carries all of cluster A's
  distinguishing loci plus a fraction `f` of cluster B's, and B gets called on genuine reads at
  the sample strain's full depth but across only `f` of its panel. These appeared at real
  abundance (1.3–34%, e.g. `S_epidermidis__GCF_006094375.1` in 11 of the mock runs), and
  identically in the 120- and 164-species panels, so neither an abundance threshold nor a
  smaller panel touches them.

  Nor can any coverage floor: measured on a shadow scenario the spurious cluster showed breadth
  0.350, while a genuinely present cluster at 0.4x showed 0.392 — indistinguishable. The two
  differ in **depth**, 7.68x versus 0.44x. Under Poisson sampling a real cluster at depth λ must
  show breadth `1 − e^(−λ)`, so the ratio is ~1 for a real cluster at any depth and ~`f` for a
  shadow. Swept synthetically: genuine clusters scored **0.949–1.018** across 0.3x–20x, shadows
  **0.200/0.300/0.495/0.691/0.897** at f = 0.2/0.3/0.5/0.7/0.9.

  Removing a shadow also repairs the true strain's abundance, which the shadow was taking a
  share of — 0.717 → 1.000 on the test scenario. So this improves composition as well as
  precision.

  Note that StrainScan's iterative "remove explained markers" step cannot be ported directly:
  cluster-unique marker sets are disjoint by construction here, so removing one cluster's
  markers leaves another's panel untouched and the iteration is a no-op. StrainScan gets its
  discrimination from scoring over the *full* k-mer set, where co-present clusters do share
  markers. The consistency test reaches the same end from the sparse-marker side, and needs no
  ordering between clusters.

  Two honest limits: the test is inert below ~0.5x depth, where `coverage ≈ 1 − e^(−depth)`
  holds for anything (so a rare strain is never penalized, but a low-depth shadow is not caught
  either); and it degrades gracefully as `f` → 1, which is correct, since a strain carrying all
  of B's distinguishing loci really is evidence for B. Set `--min-consistency 0` to disable.

## [Unreleased] — separate community members from trace contamination

### Added — `--min-global-abundance`
- Drops calls below a given share of the cross-species composition, and demotes a species whose
  clusters were all trace-level to *detected, not strain-resolvable* rather than deleting it.

  This is the **only** filter that can remove a spurious species. `--min-abundance` is
  within-species, so a species contributing a single cluster always sits at 1.0 there and can
  never be filtered by it, however little DNA it represents.

  It exists because a gate sweep on a 28-database / MSA-1005 WMS sample moved almost nothing:
  precision 0.333 at default, 0.353 under `--fixed-gate`, `--no-adaptive-singleton` or
  `--min-species-marker-frac 0.05`, and `--no-adaptive-floor` alone changed *nothing at all*.
  The 12 false-positive species turned out to be exactly the roster of a **different** mock in
  the same series (ATCC MSA-1002, *E. coli* among them), each hitting its own species-specific
  markers at ~1e-4 of the composition — three orders of magnitude below the true members at
  ~16%. That is the signature of trace cross-contamination between multiplexed libraries
  (index hopping typically leaves 0.01–0.1% of one library in another), not of a marker or
  gating error: the reads really are present and really do belong to those species.

  So no evidence gate can fix it, and tightening one only trades recall away. The gates answer
  "is there evidence"; whether an organism is a community member is a question about *how much*,
  and needs an abundance scale. On a defined mock the two populations are separated by ~3 orders
  of magnitude, so any cut inside the gap is stable — 0.001 to 0.01 is a reasonable starting
  range. Default is 0, so nothing is filtered unless asked.

## [Unreleased] — separate the two adaptive gates

### Fixed — large panels lost specificity
- **`--no-adaptive-singleton` and `--no-adaptive-floor` split what `--fixed-gate` used to
  bundle.** Real-data testing on 164-species panels showed the recall fix came with a
  precision collapse on multi-enzyme WMS samples (MSA1005/1007 precision 1.0 → 0.25–0.30;
  MSA1003 1.0 → 0.59), degrading with panel size (a 120-species panel held 0.77–0.91).

  The cause is a design flaw in the adaptive floor: it scales the species-marker floor by
  `1 − e^(−λ)` where λ is estimated from **that species' own** markers. An *absent* species has
  the lowest λ of all, so it receives the **largest** relaxation — the gate is loosest for
  exactly the species that should be rejected. An absent species then needs only
  `MIN_FLOOR_FRACTION × 200 = 50` species-specific markers hit **once** (the singleton
  admission also applies at λ≈0), where the old gate demanded 200 hit **twice**. On a large
  panel of dense multi-enzyme markers that is easy to reach by chance, and the more candidate
  species there are, the more clear it spuriously.

  The two relaxations act at different layers, which is why separating them matters:
  measured on a 0.3× single-strain sample, `--no-adaptive-singleton` loses the call while
  `--no-adaptive-floor` keeps it. Singleton admission is what buys low-depth **sensitivity**
  (it is what took MSA1002 99.9% host from 11/20 to 19/20 strain-resolved); the floor
  relaxation only moves the Layer-1 species gate, which is where the false positives enter.
  So `--no-adaptive-floor` is the setting to try first on a large panel: it should cut the
  false positives while keeping the recall gain.

  `--fixed-gate` still turns both off, for backwards compatibility.

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
