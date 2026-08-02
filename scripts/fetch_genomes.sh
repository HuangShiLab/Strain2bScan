#!/usr/bin/env bash
# Fetch a reference panel from the NCBI datasets REST API, one FASTA per accession.
#
#   fetch_genomes.sh <accession_list.tsv> <output_dir>      [JOBS=5]
#
# The list is the benchmark panel file: a header line, then one accession in column 1 per
# row. Genomes land as <output_dir>/<accession>.fna.gz, which is the naming `cluster`
# turns back into the genome name, so ground truth stated per accession joins directly
# against <out>.members.tsv.
#
# Every request bypasses any configured HTTP(S) proxy: local TLS-terminating proxies
# commonly break the handshake to *.ncbi.nlm.nih.gov, and the failure looks like a hang
# rather than an error. Resumable — an accession whose .fna.gz already exists is skipped,
# so re-running after an interruption fetches only what is missing.
#
# Keep JOBS low. NCBI throttles bursts by stalling the transfer rather than refusing it,
# and two concurrent runs of this script over the same list will race on the temp
# directories and corrupt each other's downloads.
set -uo pipefail

LIST="${1:?accession list (TSV, accession in column 1, one header line)}"
OUT="${2:?output dir}"
JOBS="${JOBS:-5}"

mkdir -p "$OUT" "$OUT/.tmp"

fetch_one() {
  local acc="$1" out="$2"
  local dest="$out/$acc.fna.gz"
  [[ -s "$dest" ]] && return 0

  local tmp="$out/.tmp/$acc"
  rm -rf "$tmp"; mkdir -p "$tmp"

  local url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
  local try
  for try in 1 2 3 4 5 6; do
    # NCBI throttles bursts by stalling the transfer rather than returning an error, so a
    # plain --max-time sits on a dead socket for its full budget. Abort instead as soon as
    # throughput stays under 5 kB/s for 20 s, and let the retry loop back off.
    if env -u HTTP_PROXY -u HTTPS_PROXY -u http_proxy -u https_proxy \
         /usr/bin/curl -sS -L --noproxy '*' --max-time 120 \
           --speed-limit 5000 --speed-time 20 -o "$tmp/g.zip" "$url" \
       && unzip -qq -o "$tmp/g.zip" -d "$tmp/x" 2>/dev/null; then
      local fna
      fna=$(find "$tmp/x" -name '*_genomic.fna' -print -quit)
      if [[ -n "$fna" && -s "$fna" ]]; then
        gzip -c "$fna" > "$dest.part" && mv "$dest.part" "$dest"
        rm -rf "$tmp"
        return 0
      fi
    fi
    sleep $((try * 2))
  done
  rm -rf "$tmp"
  echo "$acc" >> "$out/.failed"
  return 1
}
export -f fetch_one

rm -f "$OUT/.failed"
tail -n +2 "$LIST" | cut -f1 | sort -u \
  | xargs -P "$JOBS" -I{} bash -c 'fetch_one "$@"' _ {} "$OUT"

n=$(ls -1 "$OUT"/*.fna.gz 2>/dev/null | wc -l | tr -d ' ')
echo "done: $n genomes in $OUT"
[[ -s "$OUT/.failed" ]] && { echo "FAILED ($(wc -l < "$OUT/.failed" | tr -d ' ')):"; cat "$OUT/.failed"; }
rmdir "$OUT/.tmp" 2>/dev/null
exit 0
