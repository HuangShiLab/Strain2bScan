# Testing this build on HPC

Branch `strainscan-port`, at `abcaeb1`. This is what changed since the last HPC run, what is
already measured, and what still needs your machines.

Read §1 first — it says which of your existing results are still valid, so you do not re-run
things that cannot have changed.

---

## 1. What changed, and what it does to your existing numbers

| commit | change | effect on results |
|---|---|---|
| `9c5a63a` | truncated/corrupt databases are rejected instead of silently loading | none on intact databases; see §4 |
| `63c4f50` | new `--marker-source kmer` (FracMinHash sketch) | none — opt-in |
| `5aa360f` | new `batch` subcommand | none — opt-in |
| `6329f9e` | a digestion locus shared by several enzymes is counted once | **multi-enzyme databases need rebuilding**; single-enzyme unchanged |
| `84c6b8a` | k-mer windows with an N in the leading k−1 bases | k-mer mode only |
| `1ce0a79` | a `--reads` file of unknown format is rejected | none unless a file was silently misread before |
| `abcaeb1` | k-mer flags added to the usage text | none |

**`6329f9e` is the only one that moves a number, and it moves less than you would expect.**
Ten of the sixteen enzymes emit 27 bp tags and their recognition sites genuinely coincide, so a
per-enzyme scan gave one locus a copy number of 2–3 and the single-copy filter dropped it.
Measured on a 19-cluster *E. coli* panel with the 14-enzyme set:

| | before | after | truth |
|---|---|---|---|
| database markers | 166,224 | 176,734 | — |
| 70/30 mixture, abundance | 72.13 / 27.87 % | 72.22 / 27.78 % | 70 / 30 |
| 90/10 mixture, abundance | 90.95 / 9.05 % | 90.95 / 9.05 % | 90 / 10 |
| coverage | 98.46 / 93.74 % | 98.47 / 93.74 % | — |
| `support` | 4,982 / 8,375 | 5,327 / 8,907 | — |

So the panel recovers ~6% of markers that were being wrongly discarded, `support` rises with it,
and abundance moves by at most 0.09 percentage points. **Single-enzyme (BcgI) databases are
byte-identical** — 25,425 markers before and after.

Practical consequence: rebuild the multi-enzyme simulation databases when convenient, but do not
re-run abundance, recall or Bray–Curtis figures on account of this. Anything reported as
`support` or a marker count does change by ~6%.

---

## 2. Accuracy: what is measured, and what is not

**There is no full accuracy benchmark of this build.** The numbers in the paper repo were
produced by earlier versions, and the saliva set is separately invalid (§5). What follows is a
local synthetic benchmark I ran on this build — useful as a floor and as something your HPC run
can be compared against, not as a paper result.

28 draft genomes (2.4 Mb, SPAdes contigs) → 25–26 clusters. Six samples, four strains each,
log-uniform abundances spanning 3–81 %, 150 bp shotgun reads, ~1 error per 3 reads. Scored at
cluster level against known truth.

| community depth | source | recall | precision | Bray–Curtis | abundance of missed strains |
|---|---|---|---|---|---|
| 30× | BcgI | 1.000 | 1.000 | 0.006 | — |
| 30× | 14 enzymes | 1.000 | 1.000 | 0.022 | — |
| 30× | k-mer *S*=100 | 1.000 | 1.000 | 0.019 | — |
| 5× | BcgI | 1.000 | 1.000 | 0.014 | — |
| 5× | 14 enzymes | 1.000 | 1.000 | 0.023 | — |
| 5× | k-mer *S*=100 | 1.000 | 1.000 | 0.017 | — |
| 2× | BcgI | 0.818 | 1.000 | 0.036 | 3.0 – 6.0 % |
| 2× | 14 enzymes | 0.826 | 1.000 | 0.050 | 3.0 – 6.0 % |
| 2× | k-mer *S*=100 | 0.826 | 1.000 | 0.051 | 3.0 – 6.0 % |
| 1× | BcgI | 0.636 | 1.000 | 0.100 | 3.0 – 8.0 % |
| 1× | 14 enzymes | 0.565 | 1.000 | 0.151 | 3.0 – 7.4 % |
| 1× | k-mer *S*=100 | 0.565 | 1.000 | 0.141 | 3.0 – 7.4 % |

Reading it honestly:

- **Precision is 1.000 at every depth in this range** — nothing false was called, at any depth,
  under any marker source. In a separate low-depth run (community 5× and below, i.e. per-strain
  coverage under ~0.2×) precision *does* break down, to 0.33–0.83. So the no-false-positive
  behaviour holds while there is usable signal and should not be quoted as unconditional.
- **Recall is complete down to 5× community coverage** and degrades below it. Every miss is a
  strain under ~0.2× absolute coverage; nothing well-covered was ever missed.
- **The three marker sources are equivalent here.** The differences (e.g. recall 0.636 vs 0.565
  at 1×) are one strain out of 22 on six samples — noise, not a result. Do not read a ranking
  into this table.
- **This panel is easy.** 28 divergent genomes cluster into near-singletons, so any marker source
  looks good. It does not test the hard case, which is many genomes of *one* species.

Reproduce or extend it:

```bash
python3 scripts/sim_bench.py --genomes <genome_dir> --out <workdir> --bin target/release/strain2bscan
```

**The run worth doing on HPC is the same script pointed at a dense conspecific panel** — the 545
*C. acnes* references, or the 112 *P. copri* set. That is where the numbers above will not hold,
and where a real accuracy statement has to come from.

---

## 3. The k-mer marker source

`--marker-source kmer` keeps every canonical *k*-mer whose hash falls below `u64::MAX / scale`
(FracMinHash). It needs assemblies, so it is for shotgun analysis, not native 2bRAD.

Verified on *E. coli* K-12: density tracks 1/scale to within sampling noise (ratios 0.995–1.008
across scales 10–5,000), FNV-1a collides on none of the 989,962 distinct canonical 31-mers of a
1 Mb window, and the apparent shortfall against genome length is exactly the single-copy filter
removing 30,274 multi-copy 31-mers.

What it is for is **density the enzyme panel cannot reach**. On a 28-genome panel, every cluster
clears the support floor on its own markers but the Cluster Search Tree's internal nodes do not:

| source | markers/cluster (median) | node group-specific (median) | tree usable? |
|---|---|---|---|
| BcgI | 650 | 2 | no |
| `recommended` (14) | 11,711 | 32 | marginal |
| `all` (16 — the ceiling) | 16,008 | 45 | marginal |
| k-mer *S*=100 | 18,991 | 40 | marginal |
| k-mer *S*=30 | 63,554 | **134** | **yes** |
| k-mer *S*=10 | 189,986 | 372 | yes |

At matched density the two sources agree (16,008 → 45 against 18,991 → 40), so this is density,
not restriction chemistry. But `all` is the whole enzyme table and its node median is still below
the threshold to descend on. `S`=30 clears it, for ~3.5× the database size (76 MB vs 22 MB on 28
genomes; extrapolating, ~1.5 GB on 545 genomes).

**This is the likely explanation for the CST attempt failing before.** The starvation is in the
tree's *internal nodes*, not in the clusters — check with `diagnose-tree` before assuming
otherwise:

```bash
strain2bscan diagnose-tree --genomes <panel> --enzyme BcgI
strain2bscan diagnose-tree --genomes <panel> --marker-source kmer --sketch-scale 30
```

---

## 4. Check your existing databases for truncation first

Databases written before `9c5a63a` have a three-field header and can be silently truncated — a
job killed mid-write, a full quota, an interrupted transfer. A truncated database loses a
contiguous *tail* of clusters with no error, which looks like strains inexplicably going missing.

This costs seconds and needs no binary:

```bash
for f in /path/to/dbs/*.tsv; do [ "$(tail -c1 "$f" | od -An -tx1 | tr -d ' ')" = "0a" ] || echo "TRUNCATED: $f"; done
```

A complete database ends with a newline. Rebuild anything this flags, and rebuild in general to
get the marker-count check — the new header carries per-cluster marker counts and `load` verifies
them, so a truncated database becomes an error instead of a silent partial load.

Note the remaining gap: for an old three-field header, a cut landing inside the *last* row still
leaves the cluster count matching, so the binary accepts it. The shell check above catches that
case and the binary does not; there is an unpushed fix for it on the branch
`backup/local-truncation-fixes` if you want it.

---

## 5. Paper repo: the saliva analysis must be re-run

Separate from the engine, and more serious. `scripts/profile_saliva.py` parsed multi-profile's
output positionally and stopped at index 4, so `sample_fraction` was never captured; every
downstream community matrix used index 7 — `support`, a marker count, not an abundance.

Fixed in `Strain2bScan-paper@aef3cf8`: everything reads by column name now, cross-species
matrices use `sample_fraction`, per-species matrices use `within_abund`, and the scripts abort
rather than substitute a column. **They will refuse to run against the committed table**, which
is intended — it predates the fix.

```bash
# on the machine with the saliva reads
python3 scripts/profile_saliva.py      # regenerates results/saliva_strain_long.tsv
python3 scripts/analyze_saliva.py      # Fig 7 top row + PERMANOVA
```

The claim at stake: on `support`, strain-level R² 0.8328 > species-level 0.8222. On
`within_abund` the order reverses (0.6886 vs 0.7101). A 299-permutation null rules out feature
count as the explanation (null medians 0.2261 vs 0.2289 for 229 and 18 features), so the column
choice decides it. **The correct comparison on `sample_fraction` has not been run**, and no
version of this claim should go in the manuscript until it has.

See `Strain2bScan-paper/docs/rerun_status.md` for the full list of what does and does not need
re-running.

---

## 6. Suggested order

1. §4 truncation check on the existing databases — seconds, and it tells you whether the strain
   losses you saw were this.
2. §5 saliva re-run — the only thing currently blocking a manuscript number.
3. §2 benchmark on a dense conspecific panel — the real accuracy statement.
4. §3 `diagnose-tree` on the same panel, enzyme vs sketch — settles whether the CST route is
   worth reviving.
5. Rebuild multi-enzyme databases at leisure.
