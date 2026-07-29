#!/usr/bin/env python3
"""Sweep --min-global-abundance over EXISTING multi-profile .pred files.

The threshold is pure post-processing: `global_abundance` is already in the output, so a
sweep needs no re-run of the profiler. Reports species-level TP/FP/FN at each cut so the
operating point can be read off a plateau rather than guessed.

  sweep_global_abundance.py --pred run/*.pred --truth truth.tsv [--thresholds 0,1e-4,...]

truth.tsv: one expected species name per line (or "<sample>\\t<species>" for per-sample truth).
"""
import argparse, glob, os, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--pred", nargs="+", required=True)
ap.add_argument("--truth", required=True)
ap.add_argument("--thresholds", default="0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2")
ap.add_argument("--col", default="global_abundance")
a = ap.parse_args()

# truth: global list, or per-sample if two columns
truth_global, truth_by_sample = set(), defaultdict(set)
for line in open(a.truth):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    f = line.split("\t")
    (truth_by_sample[f[0]].add(f[1]) if len(f) >= 2 else truth_global.add(f[0]))

files = [p for pat in a.pred for p in glob.glob(pat)]
if not files:
    sys.exit("no .pred files matched")

def load(path):
    """-> {species: max global_abundance over its clusters}"""
    best = {}
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            si, ci = header.index("species"), header.index(a.col)
        except ValueError:
            sys.exit(f"{path}: need 'species' and '{a.col}' columns; got {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(si, ci):
                continue
            sp = f[si].strip()
            try:
                v = float(f[ci])
            except ValueError:
                continue
            best[sp] = max(best.get(sp, 0.0), v)
    return best

thresholds = [float(t) for t in a.thresholds.split(",")]
print(f"{'threshold':>10} {'TP':>5} {'FP':>5} {'FN':>5} {'precision':>10} {'recall':>8} {'F1':>7}")
for t in thresholds:
    TP = FP = FN = 0
    for path in files:
        sample = os.path.basename(path).split(".")[0]
        truth = truth_by_sample.get(sample, truth_global)
        called = {sp for sp, v in load(path).items() if v >= t}
        TP += len(called & truth); FP += len(called - truth); FN += len(truth - called)
    p = TP / (TP + FP) if TP + FP else 0.0
    r = TP / (TP + FN) if TP + FN else 0.0
    f1 = 2 * p * r / (p + r) if p + r else 0.0
    print(f"{t:>10.5g} {TP:>5} {FP:>5} {FN:>5} {p:>10.3f} {r:>8.3f} {f1:>7.3f}")
