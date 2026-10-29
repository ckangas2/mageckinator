#!/usr/bin/env python3
"""
Samplesheet + comparisons -> a validated design, written before any counting.

A pooled CRISPR screen is a subtraction between two cell populations. Which
populations, and in which direction, is the entire experiment, and MAGeCK will
run on whatever two groups it is handed without complaint. The failures this
catches all produce a gene_summary.txt that looks exactly like a real result:

  - a comparison names a condition that is not in the samplesheet: MAGeCK is
    given an empty or partial group and still writes a ranking
  - treatment and control are the same condition: every gene sits at p ~ 1,
    which reads as "a clean screen with no hits" rather than "no contrast"
  - one replicate per arm: MAGeCK estimates variance from the sgRNA mean-
    variance trend and reports confident p-values from n = 1
  - the selection direction is left implicit: in a viability/killing screen
    the hits enrich (positive selection); in a dropout screen they deplete
    (negative). Reading the wrong column of gene_summary.txt inverts the
    biology, and nothing in the file says which column you meant.

The validated design (comparisons.validated.tsv) is written once, here, and
every downstream step reads the sample lists and direction from it rather than
re-deriving them, so the comparison that runs is the comparison this report
describes.
"""
import argparse
import csv
import sys
from pathlib import Path

SELECTIONS = ("positive", "negative", "both")


def read_csv(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path} has a header but no rows")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--min-replicates", type=int, default=2,
                    help="minimum libraries per arm (plasmid baselines are "
                         "often n=1; set --min-control-replicates for those)")
    ap.add_argument("--min-control-replicates", type=int, default=None,
                    help="override the floor for the control arm only")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    min_ctrl = args.min_control_replicates
    if min_ctrl is None:
        min_ctrl = args.min_replicates

    problems = []

    samples = read_csv(args.samplesheet)
    required = {"sample", "fastq_1", "condition", "replicate"}
    missing_cols = required - set(samples[0])
    if missing_cols:
        raise SystemExit(f"samplesheet is missing columns: {sorted(missing_cols)}")

    by_condition = {}
    seen = set()
    for r in samples:
        s = r["sample"].strip()
        if s in seen:
            problems.append(f"duplicate sample name {s!r} in the samplesheet")
        seen.add(s)
        by_condition.setdefault(r["condition"].strip(), []).append(s)

    comparisons = read_csv(args.comparisons)
    required_c = {"id", "control", "treatment", "selection"}
    missing_c = required_c - set(comparisons[0])
    if missing_c:
        raise SystemExit(f"comparisons is missing columns: {sorted(missing_c)}")

    validated = []
    ids = set()
    for c in comparisons:
        cid = c["id"].strip()
        ctrl = c["control"].strip()
        trt = c["treatment"].strip()
        sel = c["selection"].strip().lower()

        if cid in ids:
            problems.append(f"duplicate comparison id {cid!r}")
        ids.add(cid)

        if sel not in SELECTIONS:
            problems.append(f"{cid}: selection {sel!r} is not one of {SELECTIONS}")

        if ctrl == trt:
            problems.append(f"{cid}: control and treatment are both {ctrl!r}; "
                            f"a comparison must contrast two different conditions")

        ctrl_samples = by_condition.get(ctrl, [])
        trt_samples = by_condition.get(trt, [])
        if not ctrl_samples:
            problems.append(f"{cid}: control condition {ctrl!r} has no samples "
                            f"in the samplesheet (have: {sorted(by_condition)})")
        if not trt_samples:
            problems.append(f"{cid}: treatment condition {trt!r} has no samples "
                            f"in the samplesheet (have: {sorted(by_condition)})")
        if ctrl_samples and len(ctrl_samples) < min_ctrl:
            problems.append(f"{cid}: control arm {ctrl!r} has {len(ctrl_samples)} "
                            f"library(ies); need {min_ctrl} "
                            f"(--min-control-replicates)")
        if trt_samples and len(trt_samples) < args.min_replicates:
            problems.append(f"{cid}: treatment arm {trt!r} has {len(trt_samples)} "
                            f"library(ies); need {args.min_replicates} "
                            f"(--min-replicates)")

        validated.append({
            "id": cid, "control": ctrl, "treatment": trt, "selection": sel,
            "control_samples": ",".join(ctrl_samples),
            "treatment_samples": ",".join(trt_samples),
            "n_control": len(ctrl_samples), "n_treatment": len(trt_samples),
        })

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1

    cols = ["id", "control", "treatment", "selection", "control_samples",
            "treatment_samples", "n_control", "n_treatment"]
    with open(out_dir / "comparisons.validated.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(validated)

    with open(out_dir / "design_report.txt", "w") as fh:
        fh.write("mageckinator design report\n")
        fh.write(f"samplesheet : {args.samplesheet}\n")
        fh.write(f"comparisons : {args.comparisons}\n\n")
        fh.write("conditions and their libraries:\n")
        for cond, ss in sorted(by_condition.items()):
            fh.write(f"  {cond:<16} n={len(ss)}  {', '.join(ss)}\n")
        fh.write("\ncomparisons:\n")
        for v in validated:
            arrow = {"positive": "enriched in", "negative": "depleted in",
                     "both": "enriched/depleted in"}[v["selection"]]
            fh.write(f"  {v['id']}: hits are guides {arrow} {v['treatment']} "
                     f"vs {v['control']} "
                     f"(n={v['n_treatment']} vs {v['n_control']})\n")

    for v in validated:
        print(f"  {v['id']:<28} {v['treatment']} vs {v['control']}  "
              f"selection={v['selection']}  "
              f"n={v['n_treatment']} vs {v['n_control']}  ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
