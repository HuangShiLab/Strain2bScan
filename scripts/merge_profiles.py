#!/usr/bin/env python3
"""Merge per-sample strain-level .pred profiles into one table.

Accepts both output layouts:
  * multi-profile --out  (#species cluster abundance coverage support depth
                          global_abundance sample_fraction n_markers)
  * profile --out        (#cluster abundance coverage support depth n_markers)

Default output is a long (tidy) TSV with one row per sample x strain call.
--matrix pivots to a wide strain x sample matrix on a chosen value column
(default: sample_fraction, the only column whose denominator is fixed by the
sequencing rather than by what the run happened to resolve — see the
multi-profile header comment in main.rs; falls back to abundance for
single-DB preds).

Sample names default to the file stem (WMS_MSA1002_0_100ng_1.pred ->
WMS_MSA1002_0_100ng_1). Ambiguous calls (Species__C0|C1) are kept verbatim.
Stdlib only.

Usage:
  merge_profiles.py --out merged.tsv *.pred
  merge_profiles.py --out matrix.tsv --matrix [--value sample_fraction] *.pred
"""
import argparse
import sys
from pathlib import Path

LONG_COLUMNS = [
    "sample", "species", "cluster", "strain_id",
    "abundance", "coverage", "support", "depth",
    "global_abundance", "sample_fraction", "n_markers",
]


def sample_name(path):
    name = Path(path).name
    for suffix in (".pred", ".tsv", ".csv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_pred(path):
    """Return (rows, is_multi). rows: list of dicts keyed by header names."""
    with open(path) as f:
        header = None
        rows = []
        for ln in f:
            if not ln.strip():
                continue
            if ln.startswith("#"):
                if header is None:
                    header = ln.rstrip("\n").lstrip("#").split("\t")
                continue
            if header is None:
                raise ValueError(f"{path}: data row before any #header")
            parts = ln.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    if header is None:
        raise ValueError(f"{path}: no header row found")
    is_multi = header[0] == "species"
    if not is_multi and header[0] != "cluster":
        raise ValueError(f"{path}: unrecognized pred header: {header}")
    return rows, is_multi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("preds", nargs="+", help="per-sample .pred files")
    ap.add_argument("--out", "-o", required=True, help="output TSV path")
    ap.add_argument("--matrix", action="store_true",
                    help="write a wide strain x sample matrix instead of a long table")
    ap.add_argument("--value", default=None,
                    help="column to pivot on with --matrix "
                         "(default: sample_fraction, or abundance if absent)")
    args = ap.parse_args()

    long_rows = []
    for path in args.preds:
        sample = sample_name(path)
        rows, is_multi = read_pred(path)
        for r in rows:
            species = r.get("species", "")
            cluster = r.get("cluster", "")
            strain_id = f"{species}__{cluster}" if species else cluster
            long_rows.append({
                "sample": sample,
                "species": species,
                "cluster": cluster,
                "strain_id": strain_id,
                "abundance": r.get("abundance", ""),
                "coverage": r.get("coverage", ""),
                "support": r.get("support", ""),
                "depth": r.get("depth", ""),
                "global_abundance": r.get("global_abundance", ""),
                "sample_fraction": r.get("sample_fraction", ""),
                "n_markers": r.get("n_markers", ""),
            })

    if not args.matrix:
        with open(args.out, "w") as w:
            w.write("\t".join(LONG_COLUMNS) + "\n")
            for r in long_rows:
                w.write("\t".join(r[c] for c in LONG_COLUMNS) + "\n")
        print(f"{len(args.preds)} samples, {len(long_rows)} rows -> {args.out}",
              file=sys.stderr)
        return

    value = args.value
    if value is None:
        value = ("sample_fraction"
                 if any(r["sample_fraction"] for r in long_rows) else "abundance")
    samples = [sample_name(p) for p in args.preds]
    strains = sorted({r["strain_id"] for r in long_rows})
    cell = {}
    for r in long_rows:
        key = (r["strain_id"], r["sample"])
        try:
            cell[key] = cell.get(key, 0.0) + float(r[value] or 0.0)
        except ValueError:
            raise SystemExit(f"--value {value}: non-numeric entry in sample "
                             f"{r['sample']} row {r['strain_id']}")
    with open(args.out, "w") as w:
        w.write("strain_id\t" + "\t".join(samples) + "\n")
        for s in strains:
            w.write(s + "\t" + "\t".join(f"{cell.get((s, smp), 0.0):.6f}"
                                         for smp in samples) + "\n")
    print(f"{len(samples)} samples x {len(strains)} strains ({value}) -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
