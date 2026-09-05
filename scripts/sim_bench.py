#!/usr/bin/env python3
"""Controlled synthetic benchmark with known truth: recall, precision and abundance error
across marker sources and sequencing depths.

    python3 scripts/sim_bench.py --genomes <dir> --out <workdir> [--bin target/release/strain2bscan]

Simulates shotgun samples from a genome directory at known strain proportions, builds a cluster
database per marker source, profiles every sample at four depths, and scores the calls against
truth at CLUSTER level (two genomes in one cluster contribute one truth entry, since that is the
resolution the tool reports at).

What this does and does not tell you. It measures the engine end to end on data whose answer is
known, which is what you want after a change. It does NOT tell you how the tool behaves on a
dense conspecific panel: a directory of divergent genomes clusters into near-singletons that are
easy to tell apart, and every marker source will look good on it. Point it at a panel of many
genomes of ONE species to see the hard case.
"""
import argparse, collections, gzip, json, os, random, subprocess, sys

RC = str.maketrans("ACGT", "TGCA")
READ_LEN = 150
STRIDES = [1, 6, 15, 30]   # full, then ~1/6, ~1/15, ~1/30 of the reads


def depth_labels(total_x):
    """Subsample strides paired with the community-wide coverage each one leaves."""
    return [(k, f"{total_x / k:g}x") for k in STRIDES]


def load_fasta(path):
    return "".join(l.strip().upper() for l in open(path) if not l.startswith(">"))


def simulate(genomes, out, n_samples, strains_per_sample, total_x, err_rate, seed):
    random.seed(seed)
    names = sorted(f.rsplit(".", 1)[0] for f in os.listdir(genomes)
                   if f.endswith((".fa", ".fna", ".fasta")))
    if len(names) < strains_per_sample:
        sys.exit(f"need >= {strains_per_sample} genomes, found {len(names)}")
    ext = {f.rsplit(".", 1)[0]: f for f in os.listdir(genomes)}
    seqs = {n: load_fasta(f"{genomes}/{ext[n]}") for n in names}
    os.makedirs(f"{out}/reads", exist_ok=True)
    truth = {}
    for si in range(n_samples):
        picked = random.sample(names, strains_per_sample)
        # log-uniform over ~2%..60%, so every sample has a genuinely rare member
        w = [10 ** random.uniform(-1.7, -0.22) for _ in picked]
        tot = sum(w)
        ab = {g: x / tot for g, x in zip(picked, w)}
        truth[f"s{si}"] = ab
        with gzip.open(f"{out}/reads/s{si}.fq.gz", "wt") as fo:
            ri = 0
            for g, a in ab.items():
                s = seqs[g]
                for _ in range(int(total_x * a * len(s) / READ_LEN)):
                    p = random.randrange(0, len(s) - READ_LEN)
                    r = s[p:p + READ_LEN]
                    if random.random() < 0.5:
                        r = r[::-1].translate(RC)
                    if random.random() < err_rate * READ_LEN:
                        q = random.randrange(READ_LEN)
                        r = r[:q] + random.choice("ACGT") + r[q + 1:]
                    fo.write(f"@r{ri}\n{r}\n+\n{'I' * READ_LEN}\n")
                    ri += 1
        print(f"  s{si}: {ri} reads  " +
              "  ".join(f"{g}={a:.3f}" for g, a in sorted(ab.items(), key=lambda x: -x[1])))
    json.dump(truth, open(f"{out}/truth.json", "w"), indent=1)
    return truth


def subsample(out, n_samples, depths):
    os.makedirs(f"{out}/reads_sub", exist_ok=True)
    for stride, label in depths:
        if stride == 1:
            continue
        for i in range(n_samples):
            with gzip.open(f"{out}/reads/s{i}.fq.gz", "rt") as fi, \
                 gzip.open(f"{out}/reads_sub/{label}_s{i}.fq.gz", "wt") as fo:
                k = 0
                while True:
                    rec = [fi.readline() for _ in range(4)]
                    if not rec[0]:
                        break
                    if k % stride == 0:
                        fo.write("".join(rec))
                    k += 1


def members(path):
    return {l.split()[0]: l.split()[1] for l in open(path) if not l.startswith("#")}


def read_pred(path):
    d = {}
    if not os.path.exists(path):
        return d
    for l in open(path):
        if l.startswith("#"):
            continue
        f = l.split("\t")
        d[f[0]] = float(f[1])
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bin", default="target/release/strain2bscan")
    ap.add_argument("--samples", type=int, default=6)
    ap.add_argument("--strains", type=int, default=4)
    ap.add_argument("--total-x", type=float, default=30)
    ap.add_argument("--error-rate", type=float, default=0.0023,
                    help="per-base substitution rate (default gives ~35%% of reads one error)")
    ap.add_argument("--seed", type=int, default=20260823)
    a = ap.parse_args()

    set14 = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
    modes = [("bcgi", ["--enzyme", "BcgI"], ["--enzyme", "BcgI"]),
             ("e14", ["--enzyme", set14], ["--enzyme", set14]),
             ("kmer", ["--marker-source", "kmer", "--sketch-scale", "100"], [])]

    os.makedirs(a.out, exist_ok=True)
    print("simulating:")
    truth = simulate(a.genomes, a.out, a.samples, a.strains, a.total_x, a.error_rate, a.seed)
    depths = depth_labels(a.total_x)
    subsample(a.out, a.samples, depths)

    for tag, build_args, _ in modes:
        subprocess.run([a.bin, "cluster", "--genomes", a.genomes, *build_args,
                        "--out", f"{a.out}/db_{tag}.tsv"],
                       check=True, stdout=subprocess.DEVNULL)
    os.makedirs(f"{a.out}/pred", exist_ok=True)
    for stride, label in depths:
        for tag, _, prof_args in modes:
            for i in range(a.samples):
                reads = (f"{a.out}/reads/s{i}.fq.gz" if stride == 1
                         else f"{a.out}/reads_sub/{label}_s{i}.fq.gz")
                subprocess.run([a.bin, "profile", "--db", f"{a.out}/db_{tag}.tsv",
                                "--reads", reads, *prof_args,
                                "--out", f"{a.out}/pred/{label}_{tag}_s{i}.tsv"],
                               check=True, stdout=subprocess.DEVNULL)

    print(f"\n{'depth':<7}{'source':<7}{'recall':>8}{'prec':>7}{'L1':>8}{'BrayCurtis':>12}"
          f"   abundance of missed strains")
    rows = []
    for _, label in depths:
        for tag, _, _ in modes:
            mm = members(f"{a.out}/db_{tag}.members.tsv")
            tp = fp = fn = 0
            l1s, bcs, missed = [], [], []
            for i in range(a.samples):
                tc = collections.defaultdict(float)
                for g, ab in truth[f"s{i}"].items():
                    tc[mm[g]] += ab
                pc = read_pred(f"{a.out}/pred/{label}_{tag}_s{i}.tsv")
                tp += len(set(tc) & set(pc))
                fn += len(set(tc) - set(pc))
                fp += len(set(pc) - set(tc))
                missed += [round(tc[c], 3) for c in set(tc) - set(pc)]
                keys = set(tc) | set(pc)
                diff = sum(abs(tc.get(c, 0) - pc.get(c, 0)) for c in keys)
                den = sum(tc.values()) + sum(pc.values())
                l1s.append(diff)
                bcs.append(diff / den if den else 1.0)
            rec = tp / (tp + fn) if tp + fn else float("nan")
            prec = tp / (tp + fp) if tp + fp else float("nan")
            print(f"{label:<7}{tag:<7}{rec:>8.3f}{prec:>7.3f}{sum(l1s)/a.samples:>8.3f}"
                  f"{sum(bcs)/a.samples:>12.3f}   {sorted(missed)[:5] if missed else '-'}")
            rows.append(dict(depth=label, source=tag, recall=rec, precision=prec,
                             l1=sum(l1s) / a.samples, bray_curtis=sum(bcs) / a.samples,
                             missed=sorted(missed)))
        print()
    json.dump(rows, open(f"{a.out}/bench_summary.json", "w"), indent=1)
    print(f"wrote {a.out}/bench_summary.json")


if __name__ == "__main__":
    main()
