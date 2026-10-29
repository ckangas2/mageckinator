#!/usr/bin/env bash
# Full path: download the GSE158298 CRISPR-screen FASTQ from SRA/ENA and write
# samplesheet.csv. Prefer scripts/fetch_counts.sh unless you want to re-count
# from reads.
#
#   bash scripts/fetch_runs.sh --list    # show how each run's title parses
#   bash scripts/fetch_runs.sh           # download screen runs + write samplesheet
#
# Accessions and URLs come from ENA's run report at run time; select_runs.py
# keeps only the pooled-screen libraries (plasmid / low-MOI / high-MOI) and
# assigns the condition from the sample title. HTTPS with aria2c where present;
# every file checked against ENA's MD5.
set -euo pipefail
cd "$(dirname "$0")/.."

PROJECT="${PROJECT:-PRJNA662661}"   # GSE158298 (verify against GEO if ENA is empty)
FASTQ_DIR="${FASTQ_DIR:-fastq}"
CONN="${CONN:-16}"
mkdir -p config "$FASTQ_DIR"

report=config/ena_runs.tsv
curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,sample_alias,sample_title,library_strategy,library_layout,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&limit=0" > "$report"
echo "ENA: $(($(wc -l < "$report") - 1)) runs in $PROJECT"

if [ "${1:-}" = "--list" ]; then
  python3 scripts/select_runs.py "$report" --list | column -t -s $'\t'
  exit 0
fi

python3 scripts/select_runs.py "$report" --fastq-dir "$FASTQ_DIR"

while IFS=$'\t' read -r u md5 f; do
  if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then
    echo "  ok (verified): $(basename "$f")"; continue; fi
  echo "  -> $(basename "$f")"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x"$CONN" -s"$CONN" -c --summary-interval=30 -d "$(dirname "$f")" -o "$(basename "$f")" "$u"
  else
    curl -fL --retry 3 -C - -o "$f" "$u"
  fi
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ] || { echo "MD5 MISMATCH: $f" >&2; exit 1; }
done < config/downloads.tsv
echo "samplesheet.csv written; FASTQs in $FASTQ_DIR/"
