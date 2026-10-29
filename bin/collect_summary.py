#!/usr/bin/env python3
"""
Run-level tables: one row per sample, one row per comparison.

Reads the count QC table and every per-comparison gene_summary and writes the
two files a reader actually opens first -- summary.tsv (how each library
sequenced) and comparisons.tsv (how many hits each comparison gave, in the
direction that comparison was run, with the top genes named). No thresholds or
gates here; this only assembles what the gated steps already decided.
"""
import argparse
import csv
import sys
from pathlib import Path


def col(row, name, default="nan"):
    v = row.get(name)
    return v if v not in (None, "") else default


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count-qc", required=True)
    ap.add_argument("--comparisons", required=True,
                    help="comparisons.validated.tsv")
    ap.add_argument("--gene-summary", nargs="*", default=[],
                    help="<id>=<path> to each gene_summary.txt")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # summary.tsv: pass the per-sample count QC straight through
    qc_rows = []
    with open(args.count_qc) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            qc_rows.append(line.rstrip("\n"))
    (out / "summary.tsv").write_text("\n".join(qc_rows) + "\n")

    selection = {}
    for r in read_tsv(args.comparisons):
        selection[r["id"]] = r["selection"]

    gs_paths = {}
    for item in args.gene_summary:
        cid, _, path = item.partition("=")
        gs_paths[cid] = path

    with open(out / "comparisons.tsv", "w") as fh:
        fh.write("id\tselection\tdirection\tn_hits_fdr\ttop_genes\n")
        for cid, path in gs_paths.items():
            rows = read_tsv(path)
            sel = selection.get(cid, "positive")
            dirs = ["pos", "neg"] if sel == "both" else \
                   (["pos"] if sel == "positive" else ["neg"])
            for d in dirs:
                hits = [r for r in rows if float(col(r, f"{d}|fdr", "1")) < args.fdr]
                hits.sort(key=lambda r: float(col(r, f"{d}|rank", "1e9")))
                top = ",".join(r["id"] for r in hits[:args.top])
                fh.write(f"{cid}\t{sel}\t{d}\t{len(hits)}\t{top}\n")
                print(f"  {cid:<28} {d}  {len(hits)} hits at FDR<{args.fdr}  "
                      f"{top[:60]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
