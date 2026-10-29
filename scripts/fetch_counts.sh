#!/usr/bin/env bash
# Fast path: fetch the deposited sgRNA read-count table for GSE158298 and write
# it as sgRNA,gene,<samples> for --count_table, skipping `mageck count`.
#
#   bash scripts/fetch_counts.sh
#
# GEO supplementary files for a series live under its acc.cgi + the
# suppl/ FTP path. The screen's count matrix is the supplementary file whose
# columns are the plasmid and MOI libraries. Column names rarely match the
# condition labels in comparisons, so after this runs, check the header and
# write samplesheet.csv (sample = each count column; condition = Plasmid /
# LowMOI / HighMOI) to match assets/gse158298_comparisons.csv.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p reference

GSE="${GSE:-GSE158298}"
stub="${GSE%???}nnn"
base="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${GSE}/suppl"

echo "GEO suppl listing for $GSE"
curl -fsSL "${base}/" -o reference/geo_suppl.html || {
  echo "could not list $base (GEO/NCBI not reachable from here)" >&2; exit 1; }

# pull candidate supplementary filenames and grab anything that looks like counts
grep -oE "${GSE}[^\"']*\.(txt|csv|tsv)(\.gz)?" reference/geo_suppl.html | sort -u \
  | while read -r name; do
      case "$name" in
        *[Cc]ount*|*sgRNA*|*CRISPR*|*readcount*)
          echo "  -> $name"
          curl -fsSL "${base}/${name}" -o "reference/${name}" ;;
      esac
    done

raw="$(ls reference/*[Cc]ount* reference/*sgRNA* 2>/dev/null | head -1 || true)"
[ -n "$raw" ] || { echo "no count-like supplementary file found; inspect reference/geo_suppl.html" >&2; exit 1; }

python3 - "$raw" reference/GSE158298_counts.tsv <<'PY'
import csv, gzip, io, re, sys
raw, out = sys.argv[1], sys.argv[2]
op = gzip.open if raw.endswith(".gz") else open
with op(raw, "rt") as fh:
    sample = fh.read()
delim = "\t" if sample.count("\t") >= sample.count(",") else ","
rows = list(csv.reader(io.StringIO(sample), delimiter=delim))
header = rows[0]
# assume col0 = sgRNA, col1 = gene when the second column is non-numeric
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sgRNA", "gene"] + header[2:])
    for r in rows[1:]:
        if len(r) >= 3:
            w.writerow([r[0], r[1]] + r[2:])
print(f"wrote {out}: {len(rows)-1} sgRNAs, columns {header[2:]}")
PY
echo
echo "Next: inspect the column names above, write samplesheet.csv mapping each"
echo "count column to a condition (Plasmid/LowMOI/HighMOI), then:"
echo "  nextflow run . -profile gse158298,docker \\"
echo "      --count_table reference/GSE158298_counts.tsv"
