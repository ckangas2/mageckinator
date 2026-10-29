#!/usr/bin/env bash
# Fetch the human GeCKOv2 sgRNA library used by GSE158298 and write it in the
# three-column form mageck count wants: sgRNA,sequence,gene.
#
#   bash scripts/fetch_library.sh
#   LIBRARY_URL=<url to a GeCKOv2 csv/xlsx> bash scripts/fetch_library.sh
#
# GeCKOv2 (Sanjana, Shalem & Zhang, Nat Methods 2014) is distributed via
# Addgene (pooled library #1000000048/#1000000049). Addgene serves the
# sequences as a zipped CSV; the exact media URL changes, so it is a parameter.
# Whatever the source columns are named, this normalises to sgRNA,sequence,gene
# and leaves the nontargeting controls labelled so validate_library finds them.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p reference

URL="${LIBRARY_URL:-}"
OUT="reference/GeCKOv2_library.csv"

if [ -z "$URL" ]; then
  cat >&2 <<'MSG'
No LIBRARY_URL set. Download the human GeCKOv2 library (Addgene
#1000000048/#1000000049, "Human GeCKO v2 library"), then either:

  LIBRARY_URL=<direct-link-to-the-csv-or-zip> bash scripts/fetch_library.sh

or drop the file at reference/GeCKOv2_raw.csv (or .xlsx) and re-run. The file
is expected to carry, per sgRNA, a guide id, a 20 nt spacer sequence, and a
gene symbol (nontargeting controls labelled "Non-Targeting Control" or similar).
MSG
  [ -f reference/GeCKOv2_raw.csv ] || [ -f reference/GeCKOv2_raw.xlsx ] || exit 1
else
  echo "-> $URL"
  tmp="reference/GeCKOv2_download"
  curl -fL --retry 3 -o "$tmp" "$URL"
  case "$URL" in
    *.zip) ( cd reference && unzip -o "$(basename "$tmp")" ) ;;
  esac
  # first csv/xlsx we can see becomes the raw input
  f="$(ls reference/*.csv reference/*.xlsx 2>/dev/null | grep -vi validated | head -1)"
  [ -n "$f" ] && cp "$f" "reference/GeCKOv2_raw.${f##*.}"
fi

python3 - "$OUT" <<'PY'
import csv, glob, re, sys
out = sys.argv[1]
src = (glob.glob("reference/GeCKOv2_raw.csv") or glob.glob("reference/GeCKOv2_raw.xlsx"))
if not src:
    sys.exit("no reference/GeCKOv2_raw.{csv,xlsx} to normalise")
src = src[0]

def rows_from(path):
    if path.endswith(".xlsx"):
        import pandas as pd
        df = pd.read_excel(path, dtype=str).fillna("")
        return [list(df.columns)] + df.values.tolist()
    with open(path, newline="") as fh:
        return list(csv.reader(fh))

rows = rows_from(src)
header = [h.strip().lower() for h in rows[0]]

def find(*names):
    for n in names:
        for i, h in enumerate(header):
            if n in h:
                return i
    return None

i_id   = find("uid", "guide", "sgrna", "name")
i_seq  = find("seq")
i_gene = find("gene")
if None in (i_seq, i_gene):
    sys.exit(f"could not find sequence/gene columns in {header}")

dna = re.compile(r"^[ACGTacgt]{18,25}$")
n = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sgRNA", "sequence", "gene"])
    for r in rows[1:]:
        if not r or len(r) <= max(i_seq, i_gene):
            continue
        seq = r[i_seq].strip().upper()
        if not dna.match(seq):
            continue
        gene = r[i_gene].strip()
        sg = r[i_id].strip() if i_id is not None else f"{gene}_{n}"
        w.writerow([sg, seq, gene])
        n += 1
print(f"wrote {n} sgRNAs to {out}")
PY
echo "library at $OUT"
