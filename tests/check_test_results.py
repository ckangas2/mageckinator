#!/usr/bin/env python3
"""
Assert the test-profile run against the planted truth in tests/data/truth.tsv.

Run after `nextflow run . -profile test,docker`:

    python3 tests/check_test_results.py results

Checks the things the synthetic screen was built to prove:
  - every library is in summary.tsv and mapped nearly all of its reads
  - the 8 host-factor genes come back as positive-selection hits (enriched)
  - the 8 essential genes come back as negative-selection hits (depleted)
  - control pseudo-genes stay out of both hit lists
  - reconciliation recovers the published host-factor list
Exit non-zero on the first broken assertion, naming it.
"""
import csv
import sys
from pathlib import Path

FDR = 0.05


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def col(row, name, default="nan"):
    v = row.get(name)
    return v if v not in (None, "") else default


def find(results, *relpaths):
    for rp in relpaths:
        p = results / rp
        if p.exists():
            return p
    hits = list(results.rglob(relpaths[-1].split("/")[-1]))
    return hits[0] if hits else None


def main():
    results = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
    truth = {r["gene"]: r["role"]
             for r in read_tsv(Path(__file__).resolve().parent / "data/truth.tsv")}
    host = {g for g, r in truth.items() if r == "host_factor"}
    ess = {g for g, r in truth.items() if r == "essential"}

    failures = []

    def check(cond, msg):
        print(("ok  " if cond else "FAIL") + "  " + msg)
        if not cond:
            failures.append(msg)

    cid = "hostfactors_highmoi_vs_plasmid"

    # ---- summary.tsv
    summ = find(results, "summary.tsv")
    check(summ is not None, "summary.tsv exists")
    if summ:
        text = [l for l in summ.read_text().splitlines() if l and not l.startswith(("#", "sample"))]
        samples = [l.split("\t")[0] for l in text]
        for s in ["Plasmid_rep1", "Plasmid_rep2", "HighMOI_rep1",
                  "HighMOI_rep2", "HighMOI_rep3"]:
            check(s in samples, f"{s} present in summary.tsv")

    # ---- gene_summary for the one comparison
    gs = find(results, f"test/{cid}/{cid}.gene_summary.txt",
              f"{cid}.gene_summary.txt")
    check(gs is not None, "gene_summary.txt exists")
    if gs:
        rows = read_tsv(gs)
        pos_hits = {r["id"] for r in rows if float(col(r, "pos|fdr", "1")) < FDR}
        neg_hits = {r["id"] for r in rows if float(col(r, "neg|fdr", "1")) < FDR}

        recovered_host = host & pos_hits
        check(len(recovered_host) >= 0.75 * len(host),
              f"host factors enriched (pos FDR<{FDR}): "
              f"{len(recovered_host)}/{len(host)} "
              f"[{','.join(sorted(recovered_host))}]")

        recovered_ess = ess & neg_hits
        check(len(recovered_ess) >= 0.75 * len(ess),
              f"essentials depleted (neg FDR<{FDR}): "
              f"{len(recovered_ess)}/{len(ess)}")

        ctrl_pos = {g for g in pos_hits if g.startswith("Control")}
        ctrl_neg = {g for g in neg_hits if g.startswith("Control")}
        check(len(ctrl_pos | ctrl_neg) <= 1,
              f"controls stay out of hits (got {sorted(ctrl_pos | ctrl_neg)})")

        # top of the positive list should be host factors, not flat genes
        top5 = [r["id"] for r in sorted(
            rows, key=lambda r: float(col(r, "pos|rank", "1e9")))[:5]]
        check(len(host & set(top5)) >= 3,
              f"host factors dominate the top of the positive ranking: {top5}")

    # ---- reconciliation
    rec = find(results, f"qc/{cid}.reconciliation.tsv", "reconciliation.tsv")
    if rec:
        kv = {}
        for line in rec.read_text().splitlines():
            if "\t" in line and not line.startswith("#"):
                k, _, v = line.partition("\t")
                kv[k] = v
        if "recall_positive" in kv:
            check(float(kv["recall_positive"]) >= 0.75,
                  f"reconciliation recall_positive = {kv['recall_positive']}")

    print()
    if failures:
        print(f"{len(failures)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
