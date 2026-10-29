#!/usr/bin/env python3
"""
Our ranking against the published one.

A reanalysis is only worth anything if it lands near the paper it reanalyses,
and the comparison has to be fixed in advance or it becomes a story told after
the fact. This reads a MAGeCK gene_summary and a list of the genes the original
study reported, and asks two questions:

  - recall: of the published hits, how many does our ranking recover (in the
    top-N, or below our FDR)? Below a floor, we have not reproduced the screen.
  - direction: does the published hit set agree with the *selection* column we
    chose (positive = enriched, negative = depleted)? If the opposite column
    recovers the published hits better, treatment and control were almost
    certainly swapped somewhere upstream -- a mislabel that MAGeCK cannot see,
    because the file is a valid ranking either way.

Perfect agreement is not expected: a different aligner, a different MAGeCK
version, and count vs count-table starting points all move ranks around. The
floor is for catching an inverted or scrambled comparison, not for grading the
fourth decimal place.
"""
import argparse
import csv
import math
import sys


def read_gene_summary(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def col(row, name, default="nan"):
    v = row.get(name)
    return v if v not in (None, "") else default


def load_published(path):
    with open(path, newline="") as fh:
        head = fh.readline()
        delim = "\t" if "\t" in head else ","
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter=delim)
        genes, ranks = [], {}
        for r in rdr:
            g = (r.get("gene") or r.get("id") or r.get("Gene") or "").strip()
            if not g:
                continue
            genes.append(g)
            rk = r.get("rank") or r.get("Rank")
            if rk not in (None, ""):
                ranks[g] = float(rk)
    return genes, ranks


def ranked_ids(rows, direction):
    key = f"{direction}|rank"
    if key in rows[0]:
        return [r["id"] for r in sorted(rows, key=lambda r: float(col(r, key, "1e9")))]
    score = f"{direction}|score"
    return [r["id"] for r in sorted(rows, key=lambda r: float(col(r, score, "1")))]


def recall(rows, published, direction, fdr, top_n):
    pub = set(published)
    present = pub & {r["id"] for r in rows}
    if not present:
        return 0.0, 0, set()
    top = set(ranked_ids(rows, direction)[:top_n])
    fdr_hits = {r["id"] for r in rows if float(col(r, f"{direction}|fdr", "1")) < fdr}
    recovered = present & (top | fdr_hits)
    return len(recovered) / len(present), len(present), recovered


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-summary", required=True)
    ap.add_argument("--published", required=True,
                    help="gene[,rank] list the original study reported")
    ap.add_argument("--selection", choices=("positive", "negative"),
                    default="positive")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--min-recall", type=float, default=0.5)
    ap.add_argument("--out", default="reconciliation.tsv")
    args = ap.parse_args()

    rows = read_gene_summary(args.gene_summary)
    if not rows:
        raise SystemExit("gene_summary has no rows")
    published, _pub_ranks = load_published(args.published)
    if not published:
        raise SystemExit("no published genes loaded")

    direction = {"positive": "pos", "negative": "neg"}[args.selection]
    opposite = "neg" if direction == "pos" else "pos"

    rec, n_present, recovered = recall(rows, published, direction, args.fdr, args.top_n)
    rec_opp, _, recovered_opp = recall(rows, published, opposite, args.fdr, args.top_n)

    problems = []
    print(f"  published hits present in our matrix : {n_present}/{len(set(published))}")
    print(f"  recovered in top-{args.top_n} or FDR<{args.fdr} ({args.selection}) : "
          f"{len(recovered)} ({rec:.1%})")
    print(f"  same, opposite direction             : "
          f"{len(recovered_opp)} ({rec_opp:.1%})")

    if n_present == 0:
        problems.append("none of the published genes are in our gene_summary; "
                        "gene ids do not match the library")
    else:
        if rec < args.min_recall:
            problems.append(f"recovered {rec:.1%} of published hits "
                            f"(min {args.min_recall:.0%}): this is not the same "
                            f"screen")
        if rec_opp > rec and rec_opp >= args.min_recall:
            problems.append(f"the opposite selection recovers the published hits "
                            f"better ({rec_opp:.1%} vs {rec:.1%}): treatment and "
                            f"control look swapped, or --selection is wrong")

    with open(args.out, "w") as fh:
        fh.write(f"# selection={args.selection} top_n={args.top_n} fdr={args.fdr}\n")
        fh.write(f"published_present\t{n_present}\n")
        fh.write(f"recall_{args.selection}\t{rec:.4f}\n")
        fh.write(f"recall_opposite\t{rec_opp:.4f}\n\n")
        fh.write("gene\trecovered\tour_rank\tour_fdr\n")
        rank_of = {g: i + 1 for i, g in enumerate(ranked_ids(rows, direction))}
        fdr_of = {r["id"]: col(r, f"{direction}|fdr", "nan") for r in rows}
        for g in published:
            if g in rank_of:
                fh.write(f"{g}\t{'yes' if g in recovered else 'no'}\t"
                         f"{rank_of[g]}\t{fdr_of.get(g, 'nan')}\n")
            else:
                fh.write(f"{g}\tabsent\tNA\tNA\n")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
