#!/usr/bin/env bash
# Profile every sample under every Layer-1/Layer-2 combination, against ONE shared cluster DB.
#
#   BASELINE=path/to/pre-port/strain2bscan \
#   PORT=target/release/strain2bscan \
#   DB=acnes.db.tsv READS=reads/ OUT=out/ \
#   run_bench.sh
#
# Reads are <READS>/<sample>.fq.gz for each name in SAMPLES; predictions land in
# <OUT>/<arm>/<sample>.tsv, which is the layout scripts/evaluate.py expects.
#
# Why one DB for both binaries. Only a binary carrying the port serializes the Cluster
# Search Tree, and `--layer1 cst` falls back to the flat path without it, so the DB must be
# built by $PORT (`cluster`, not `build`). The pre-port loader skips '#'-prefixed sections
# it does not recognise, so it reads that same file correctly and simply ignores the tree.
# That is what makes this a controlled comparison: identical clusters, identical marker
# sets, only the algorithm differs. Building a DB per binary instead would confound an
# algorithm difference with a clustering difference.
#
# The DB must be rebuilt when leaf marker sets change; a stale DB will not error, it will
# quietly profile against the old panels.
#
# BASELINE is optional. Set it to a binary built from the pre-port commit to get the
# regression control (its output must be byte-identical to $PORT's default path); leave it
# unset to run the port's arms alone.
set -euo pipefail

PORT="${PORT:-target/release/strain2bscan}"
BASELINE="${BASELINE:-}"
DB="${DB:?set DB to a cluster DB built by \$PORT}"
READS="${READS:?set READS to the directory holding <sample>.fq.gz}"
OUT="${OUT:-out}"
LOGS="${LOGS:-$OUT/logs}"
SAMPLES="${SAMPLES:-sample1 sample2 sample3 sample4 sample5}"

mkdir -p "$LOGS"

run() {          # run <arm> <binary> [extra profile flags...]
  local arm="$1" bin="$2"; shift 2
  mkdir -p "$OUT/$arm"
  for s in $SAMPLES; do
    "$bin" profile --db "$DB" --reads "$READS/$s.fq.gz" \
      --out "$OUT/$arm/$s.tsv" "$@" > "$LOGS/$arm.$s.log" 2>&1
  done
  echo "  $arm"
}

echo "profiling into $OUT:"
[ -n "$BASELINE" ] && run baseline-flat "$BASELINE"
run port-flat     "$PORT"
run port-cst      "$PORT" --layer1 cst
run port-enet     "$PORT" --layer2 enet
run port-cst-enet "$PORT" --layer1 cst --layer2 enet

# The control worth checking first: with both layers on their defaults the port must not
# have changed anything. If these differ, the comparison below is measuring a regression
# rather than the layers.
if [ -n "$BASELINE" ]; then
  echo
  for s in $SAMPLES; do
    cmp -s "$OUT/baseline-flat/$s.tsv" "$OUT/port-flat/$s.tsv" \
      || echo "CONTROL FAILED: baseline and port default differ on $s"
  done
  echo "control: baseline vs port default checked on $(echo $SAMPLES | wc -w | tr -d ' ') samples"
fi
