#!/usr/bin/env python3
"""Score strain2bscan predictions against a ground truth.

Two truth layouts and two prediction layouts, in any combination:

  --truth-format wide   MockMetagenomes4Benchmark's genome x sample matrix. Truth is per
                        *genome*, so it is remapped through the cluster membership sidecar
                        (--members) and summed within a cluster.
  --truth-format long   one file per sample, carrying its own cluster column. Single-species
                        files are keyed by cluster; multi-species files (a `species` column)
                        by species+cluster.

  --pred-format profile        `#cluster abundance ...`
  --pred-format multi-profile  `#species cluster abundance ... global_abundance ...`

`--abundance-col` must match how the truth is normalised, and nothing will warn you if it
does not: detection metrics are identical either way. On a three-cluster check where the
truth is community-scope, scoring the within-species column instead left precision, recall
and F1 at 1.000 while Bray-Curtis went from 0.010 to 0.333. Single-species sim truth sums to
1 WITHIN a species -> `abundance`; multi-species sim truth sums to 1 across the community ->
`global_abundance`; a truth stated as DNA/sequence mass rather than cell counts ->
`sample_fraction`. The defaults follow --pred-format.

A run is
scored on presence (precision/recall/F1 at a threshold, plus AUPR over the ranked
prediction) and on composition, the latter twice: L1/Bray-Curtis over the label union,
which is the end-to-end score and is dominated by detection, and BC_tp/MARE_tp over the
clusters both agree on, which isolates abundance accuracy. A change that only touches
abundance estimation -- Layer-2 -- moves the second pair and nothing else, so reporting
only the first pair makes it look inert.

AUPR is average precision over the predicted list ranked by abundance: sweep the
presence threshold down through the predictions and take sum_k (R_k - R_{k-1}) * P_k.
A truth cluster the run never emitted is unreachable at any threshold, so it caps
recall -- which is the intended behaviour, since an omitted cluster is a miss no
threshold can recover.
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict


def read_members(path):
    """genome -> cluster id."""
    m = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            g, c = line.rstrip("\n").split("\t")[:2]
            m[g] = c
    return m


def read_truth(path):
    """sample label -> {genome: abundance}, from the wide benchmark matrix."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    samples = header[1:]
    out = {s: {} for s in samples}
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        genome = r[0].strip()
        for s, v in zip(samples, r[1:]):
            v = float(v)
            if v > 0:
                out[s][genome] = v
    return out


def read_truth_long(path):
    """One sample's truth from a long-format table -> {label: abundance}.

    Two shapes, told apart by the header rather than by a flag:

      single  strain, cluster, relative_abundance, ...   -> label = cluster
      multi   species, strain, cluster, relative_abundance, ... -> label = "species\\tcluster"

    Truth is stated per *strain* and several strains can share a cluster, so abundances are
    summed into the cluster — the same remap `read_members` does for the wide format, except
    the mapping is already in the file.

    The scope differs between the two and the caller must match the prediction column to it:
    single-species truth is normalised WITHIN a species (the k strains sum to 1), multi-species
    truth is normalised across the whole community. See `--abundance-col`.
    """
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    if not rows:
        return {}, False
    header = [h.lstrip("#").strip() for h in rows[0]]
    idx = {name: i for i, name in enumerate(header)}
    for need in ("cluster", "relative_abundance"):
        if need not in idx:
            sys.exit(f"{path}: long truth needs a '{need}' column, got {header}")
    multi = "species" in idx
    out = defaultdict(float)
    for r in rows[1:]:
        if not r or not r[0].strip() or r[0].startswith("#"):
            continue
        try:
            v = float(r[idx["relative_abundance"]])
        except (ValueError, IndexError):
            continue
        if v <= 0:
            continue
        cluster = r[idx["cluster"]].strip()
        key = f"{r[idx['species']].strip()}\t{cluster}" if multi else cluster
        out[key] += v
    return dict(out), multi


def _species_key(raw):
    """Normalise `multi-profile`'s species column to the truth's species name.

    The column holds the DB file's stem, so a database at `Cutibacterium_acnes.db.tsv` reports
    `Cutibacterium_acnes.db` while the truth says `Cutibacterium_acnes`. Left alone, every
    label mismatches and the run scores zero across the board — which looks like a catastrophic
    result rather than a join failure.
    """
    return raw[:-3] if raw.endswith(".db") else raw


def read_pred(path, clade_mode="strict", pred_format="profile", abundance_col=None):
    """label -> abundance, from either subcommand's prediction file.

    The two subcommands write different layouts, and reading one as the other fails silently:

      profile        #cluster abundance coverage support depth n_markers
      multi-profile  #species cluster abundance coverage support depth global_abundance
                     sample_fraction n_markers

    Fed a multi-profile file, the single-species reader takes the SPECIES name as the cluster
    and its second column as the abundance, producing labels that match nothing. Columns are
    therefore located by header name, not by position, and `abundance_col` picks which of the
    three scopes to score — they answer different questions and the truth's normalisation
    decides which one is correct.

    --layer1 cst can resolve a strain only to a clade and report it under every leaf it
    spans, as `C1|C3`. Two ways to score that, and they answer different questions:

      strict  keep `C1|C3` as its own label. It matches no truth cluster, so an
              unresolved clade counts as one FP plus an FN for each true member. This is
              the right reading if the deliverable is a named strain.
      split   expand into its leaves at abundance/n_leaves each. Gives the clade call
              credit for the members it does contain, at the cost of inventing the
              siblings it does not. This is the right reading if a clade-level answer is
              useful downstream.

    Neither is more correct in general, so both are reported.
    """
    out = defaultdict(float)
    if not os.path.exists(path):
        return dict(out)
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    if not rows:
        return dict(out)
    header = [h.lstrip("#").strip() for h in rows[0]]
    idx = {name: i for i, name in enumerate(header)}
    col = abundance_col or ("global_abundance" if pred_format == "multi-profile" else "abundance")
    if col not in idx:
        sys.exit(f"{path}: no '{col}' column; header is {header}")
    if pred_format == "multi-profile" and "species" not in idx:
        sys.exit(f"{path}: --pred-format multi-profile but no 'species' column; header is {header}")
    for r in rows[1:]:
        if not r or not r[0].strip() or r[0].startswith("#"):
            continue
        try:
            ab = float(r[idx[col]])
        except (ValueError, IndexError):
            continue
        cluster = r[idx["cluster"]].strip()
        prefix = f"{_species_key(r[idx['species']].strip())}\t" if pred_format == "multi-profile" else ""
        leaves = cluster.split("|")
        if clade_mode == "split" and len(leaves) > 1:
            for leaf in leaves:
                out[f"{prefix}{leaf}"] += ab / len(leaves)
        else:
            out[f"{prefix}{cluster}"] += ab
    return dict(out)


def remap_truth(truth_genomes, members):
    """Sum genome-level truth into cluster-level truth. Returns (dict, unmapped list)."""
    out = defaultdict(float)
    unmapped = []
    for g, v in truth_genomes.items():
        c = members.get(g)
        if c is None:
            unmapped.append(g)
            continue
        out[c] += v
    return dict(out), unmapped


def _mean(vals):
    """Mean over the non-NaN entries — NaN marks a sample with nothing to score, which is
    not the same as a sample scoring zero and must not be averaged in as one."""
    v = [x for x in vals if x == x]
    return sum(v) / len(v) if v else float("nan")


def prf(pred, truth, thresh):
    p = {k for k, v in pred.items() if v >= thresh}
    t = set(truth)  # truth entries are all > 0 by construction
    tp, fp, fn = len(p & t), len(p - t), len(t - p)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return dict(tp=tp, fp=fp, fn=fn, precision=prec, recall=rec, f1=f1,
                fp_labels=sorted(p - t), fn_labels=sorted(t - p))


def aupr(pred, truth):
    """Average precision over predictions ranked by abundance descending."""
    t = set(truth)
    if not t:
        return float("nan")
    ranked = sorted(pred.items(), key=lambda kv: -kv[1])
    tp = 0
    prev_recall = 0.0
    ap = 0.0
    for i, (name, _) in enumerate(ranked, start=1):
        if name in t:
            tp += 1
        recall = tp / len(t)
        precision = tp / i
        ap += (recall - prev_recall) * precision
        prev_recall = recall
    return ap


def composition(pred, truth):
    """Abundance error, over the label union AND over the detected set alone.

    Both are needed, and for different questions. The union metrics (l1, bray_curtis) are
    the honest end-to-end score: a missed cluster contributes its full truth abundance as
    error, so detection failures show up here. But that also means they are dominated by
    detection, which makes them useless for judging a change that only touches abundance —
    Layer-2 never alters *which* clusters are called, so any effect it has is invisible in
    a metric that detection differences swamp.

    The TP-restricted metrics answer the other question: given the clusters both agree on,
    how well is the split estimated? Renormalising over that shared set first is what makes
    the comparison fair — otherwise a run that missed a 30% cluster is scored on a
    prediction vector that cannot sum to 1 against a truth vector that does.
    """
    labels = set(pred) | set(truth)
    l1 = num = den = 0.0
    for l in labels:
        pv, tv = pred.get(l, 0.0), truth.get(l, 0.0)
        l1 += abs(pv - tv)
        num += abs(pv - tv)
        den += pv + tv
    out = dict(l1=l1, bray_curtis=(num / den if den else 0.0))

    shared = sorted(set(pred) & set(truth))
    if not shared:
        out.update(l1_tp=float("nan"), bc_tp=float("nan"), mare_tp=float("nan"))
        return out
    ps, ts = sum(pred[l] for l in shared), sum(truth[l] for l in shared)
    l1_tp = num_tp = den_tp = mare = 0.0
    for l in shared:
        pv = pred[l] / ps if ps else 0.0
        tv = truth[l] / ts if ts else 0.0
        l1_tp += abs(pv - tv)
        num_tp += abs(pv - tv)
        den_tp += pv + tv
        if tv > 0:
            mare += abs(pv - tv) / tv
    out.update(l1_tp=l1_tp,
               bc_tp=(num_tp / den_tp if den_tp else 0.0),
               mare_tp=mare / len(shared))
    return out


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--members",
                     help="genome->cluster sidecar; required for --truth-format wide only "
                          "(long truth already carries a cluster column)")
    ap_.add_argument("--truth", required=True,
                     help="wide: the truth matrix file. long: the directory of per-sample "
                          "<sample><--truth-suffix> files")
    ap_.add_argument("--truth-format", choices=["wide", "long"], default="wide")
    ap_.add_argument("--truth-suffix", default=".truth.tsv",
                     help="long only: per-sample truth filename suffix")
    ap_.add_argument("--pred-format", choices=["profile", "multi-profile"], default="profile")
    ap_.add_argument("--abundance-col",
                     choices=["abundance", "global_abundance", "sample_fraction"],
                     help="which prediction column to score. Default: `abundance` for "
                          "profile, `global_abundance` for multi-profile. Match this to how "
                          "the truth is normalised — single-species sim truth sums to 1 "
                          "WITHIN a species (abundance), multi-species sim truth sums to 1 "
                          "across the community (global_abundance). Use sample_fraction only "
                          "if the truth is a DNA/sequence fraction rather than a cell one.")
    ap_.add_argument("--preds", required=True,
                     help="dir holding <run>/<sample>.tsv")
    ap_.add_argument("--samples", default="sample1,sample2,sample3,sample4,sample5")
    ap_.add_argument("--truth-samples", default="sample_1,sample_2,sample_3,sample_4,sample_5",
                     help="matching column names in the truth matrix")
    ap_.add_argument("--thresh", type=float, default=0.01)
    ap_.add_argument("--clade-mode", choices=["strict", "split"], default="strict")
    ap_.add_argument("--json-out")
    args = ap_.parse_args()

    samples = args.samples.split(",")
    tsamples = args.truth_samples.split(",")
    if args.truth_format == "wide":
        if not args.members:
            sys.exit("--members is required for --truth-format wide")
        members = read_members(args.members)
        truth_all = read_truth(args.truth)
    else:
        # Long truth is one file per sample and already carries the cluster assignment, so no
        # membership sidecar is needed.
        members, truth_all = None, {}
        for s_ in samples:
            path = os.path.join(args.truth, f"{s_}{args.truth_suffix}")
            if not os.path.exists(path):
                sys.exit(f"missing truth file {path}")
            truth_all[s_], is_multi = read_truth_long(path)
            if is_multi and args.pred_format != "multi-profile":
                sys.exit(f"{path} is multi-species truth (it has a `species` column) but "
                         "--pred-format is `profile`; the labels cannot join")
            if not is_multi and args.pred_format == "multi-profile":
                sys.exit(f"{path} is single-species truth but --pred-format is "
                         "`multi-profile`; the labels cannot join")
        tsamples = samples

    # A run directory is one holding at least one <sample>.tsv. Testing for that rather than
    # taking every subdirectory keeps siblings like the profiler's logs/ out of the table —
    # they would otherwise score 0 across the board and read as a catastrophically bad arm.
    def is_run(d):
        p = os.path.join(args.preds, d)
        return os.path.isdir(p) and any(
            os.path.exists(os.path.join(p, f"{s}.tsv")) for s in samples
        )

    runs = sorted(d for d in os.listdir(args.preds) if is_run(d))
    if not runs:
        sys.exit(f"no run subdirectories with <sample>.tsv under {args.preds}")

    results = {}
    for run in runs:
        per_sample = []
        for s, ts in zip(samples, tsamples):
            if args.truth_format == "wide":
                truth_c, unmapped = remap_truth(truth_all[ts], members)
            else:
                truth_c, unmapped = truth_all[ts], []
            path = os.path.join(args.preds, run, f"{s}.tsv")
            pred = read_pred(path, args.clade_mode, args.pred_format, args.abundance_col)
            raw = read_pred(path, "strict", args.pred_format, args.abundance_col)
            row = dict(sample=s,
                       n_truth_clusters=len(truth_c),
                       n_truth_genomes=len(truth_all[ts]),
                       n_pred=len(pred),
                       n_clade_calls=sum(1 for k in raw if "|" in k),
                       unmapped_truth_genomes=unmapped)
            row.update(prf(pred, truth_c, args.thresh))
            row["aupr"] = aupr(pred, truth_c)
            row.update(composition(pred, truth_c))
            per_sample.append(row)

        n = len(per_sample)
        agg = dict(
            precision=sum(r["precision"] for r in per_sample) / n,
            recall=sum(r["recall"] for r in per_sample) / n,
            f1=sum(r["f1"] for r in per_sample) / n,
            aupr=sum(r["aupr"] for r in per_sample) / n,
            l1=sum(r["l1"] for r in per_sample) / n,
            bray_curtis=sum(r["bray_curtis"] for r in per_sample) / n,
            # Abundance-only metrics average over the samples that HAVE a shared label;
            # a sample where nothing was detected has no abundance error to report.
            bc_tp=_mean([r["bc_tp"] for r in per_sample]),
            mare_tp=_mean([r["mare_tp"] for r in per_sample]),
            tp=sum(r["tp"] for r in per_sample),
            fp=sum(r["fp"] for r in per_sample),
            fn=sum(r["fn"] for r in per_sample),
        )
        # Micro-average: pool the counts rather than averaging per-sample rates, so a
        # sample with few truth clusters does not carry the same weight as a rich one.
        agg["micro_precision"] = agg["tp"] / (agg["tp"] + agg["fp"]) if agg["tp"] + agg["fp"] else 0.0
        agg["micro_recall"] = agg["tp"] / (agg["tp"] + agg["fn"]) if agg["tp"] + agg["fn"] else 0.0
        mp, mr = agg["micro_precision"], agg["micro_recall"]
        agg["micro_f1"] = 2 * mp * mr / (mp + mr) if mp + mr else 0.0
        results[run] = dict(per_sample=per_sample, aggregate=agg)

    w = 26
    col = args.abundance_col or (
        "global_abundance" if args.pred_format == "multi-profile" else "abundance")
    print(f"presence threshold = {args.thresh}   clade calls scored as: {args.clade_mode}")
    print(f"truth={args.truth_format}  pred={args.pred_format}  scoring column: {col}\n")
    hdr = (f"{'run':<{w}}{'P':>7}{'R':>7}{'F1':>7}{'AUPR':>8}"
           f"{'microP':>8}{'microR':>8}{'TP':>5}{'FP':>5}{'FN':>5}"
           f"{'L1':>8}{'BC':>7}{'BC_tp':>8}{'MARE_tp':>9}")
    print(hdr)
    print("-" * len(hdr))
    for run, r in results.items():
        a = r["aggregate"]
        print(f"{run:<{w}}{a['precision']:>7.3f}{a['recall']:>7.3f}{a['f1']:>7.3f}"
              f"{a['aupr']:>8.3f}{a['micro_precision']:>8.3f}{a['micro_recall']:>8.3f}"
              f"{a['tp']:>5}{a['fp']:>5}{a['fn']:>5}{a['l1']:>8.3f}{a['bray_curtis']:>7.3f}"
              f"{a['bc_tp']:>8.3f}{a['mare_tp']:>9.3f}")

    print("\nper-sample detail")
    for run, r in results.items():
        print(f"\n  [{run}]")
        for s in r["per_sample"]:
            print(f"    {s['sample']}  truth={s['n_truth_clusters']}c/{s['n_truth_genomes']}g"
                  f"  pred={s['n_pred']}  TP={s['tp']} FP={s['fp']} FN={s['fn']}"
                  f"  P={s['precision']:.3f} R={s['recall']:.3f} AUPR={s['aupr']:.3f}"
                  f"  BC={s['bray_curtis']:.3f}"
                  + (f"  clade_calls={s['n_clade_calls']}" if s["n_clade_calls"] else ""))
            if s["fp_labels"]:
                print(f"        FP: {','.join(s['fp_labels'])}")
            if s["fn_labels"]:
                print(f"        FN: {','.join(s['fn_labels'])}")
            if s["unmapped_truth_genomes"]:
                print(f"        !! truth genomes absent from the panel: "
                      f"{','.join(s['unmapped_truth_genomes'])}")

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(results, fh, indent=2)


if __name__ == "__main__":
    main()
