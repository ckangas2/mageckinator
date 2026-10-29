#!/usr/bin/env python3
"""
Gate on the count table before any gene-level test.

`mageck count` always writes a countsummary and a count matrix, whatever went
in. Each number below has a plausible-looking failure mode:

  - mapping rate: reads that do not match any guide are dropped silently. A
    wrong library, the wrong read, or un-trimmed adapter gives a full count
    table built from the few percent of reads that happened to match.
  - zero-count fraction / coverage: too few reads per guide and the counts are
    Poisson noise. MAGeCK still ranks them; the ranking is of the noise.
  - Gini index: a bottleneck or a PCR jackpot concentrates reads in a handful
    of guides. The library is no longer a fair census, but the matrix is full.
  - replicate structure: a swapped or mislabelled plasmid/baseline library is
    invisible in its own counts. It shows up as a library that correlates with
    the wrong condition's replicates. Survivor arms are expected to diverge
    (that is the signal), so this check is applied only to the control arms,
    where replicates should still agree.

Thresholds are deliberately conservative defaults; set them per screen in the
profile. The Gini and zero-count gates apply to control (baseline/plasmid)
arms by default, since positive-selection survivors legitimately skew.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np


def read_countsummary(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def read_counts(path):
    with open(path, newline="") as fh:
        delim = "\t"
        header = fh.readline().rstrip("\n").split(delim)
        samples = header[2:]
        mat = {s: [] for s in samples}
        sgrnas = []
        for line in fh:
            f = line.rstrip("\n").split(delim)
            sgrnas.append(f[0])
            for s, v in zip(samples, f[2:]):
                mat[s].append(float(v))
    return samples, sgrnas, {s: np.asarray(v) for s, v in mat.items()}


def spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = math.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / denom) if denom else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--countsummary", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--comparisons", required=True,
                    help="comparisons.validated.tsv from check_design")
    ap.add_argument("--min-mapped-frac", type=float, default=0.6)
    ap.add_argument("--max-zero-frac", type=float, default=0.3)
    ap.add_argument("--max-gini", type=float, default=0.35)
    ap.add_argument("--min-coverage", type=float, default=100.0,
                    help="mapped reads per sgRNA")
    ap.add_argument("--min-replicate-rho", type=float, default=0.8)
    ap.add_argument("--skew-arms", choices=("control", "all"), default="control",
                    help="which arms the Gini/zero-frac gates apply to")
    ap.add_argument("--out", default="count_qc.tsv")
    args = ap.parse_args()

    with open(args.samplesheet, newline="") as fh:
        cond = {r["sample"].strip(): r["condition"].strip()
                for r in csv.DictReader(fh)}

    control_conditions = set()
    with open(args.comparisons, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            control_conditions.add(r["control"].strip())

    summ = read_countsummary(args.countsummary)
    samples, sgrnas, mat = read_counts(args.counts)
    n_sgrnas = len(sgrnas)
    problems = []

    # per-sample QC from the countsummary
    qc = {}
    for r in summ:
        label = r.get("Label", r.get("label", "")).strip()
        reads = float(r.get("Reads", 0) or 0)
        mapped = float(r.get("Mapped", 0) or 0)
        pct = r.get("Percentage")
        mapped_frac = (float(pct) if pct not in (None, "") else
                       (mapped / reads if reads else 0.0))
        total = float(r.get("TotalsgRNAs", n_sgrnas) or n_sgrnas)
        zero = float(r.get("Zerocounts", 0) or 0)
        gini = float(r.get("GiniIndex", "nan") or "nan")
        qc[label] = dict(reads=reads, mapped=mapped, mapped_frac=mapped_frac,
                         zero_frac=(zero / total if total else 0.0),
                         coverage=(mapped / total if total else 0.0), gini=gini)

    for s in samples:
        c = cond.get(s)
        if c is None:
            problems.append(f"{s} is in the count table but not the samplesheet")
        if s not in qc:
            problems.append(f"{s} is in the count table but not the countsummary")

    def arm_in_scope(sample):
        return args.skew_arms == "all" or cond.get(sample) in control_conditions

    for s in samples:
        q = qc.get(s)
        if q is None:
            continue
        if q["mapped_frac"] < args.min_mapped_frac:
            problems.append(f"{s}: {q['mapped_frac']:.1%} of reads map to the "
                            f"library (min {args.min_mapped_frac:.0%}): wrong "
                            f"library, wrong read, or un-trimmed adapter")
        if q["coverage"] < args.min_coverage:
            problems.append(f"{s}: {q['coverage']:.0f} mapped reads per sgRNA "
                            f"(min {args.min_coverage:.0f}): under-sequenced")
        if arm_in_scope(s):
            if q["zero_frac"] > args.max_zero_frac:
                problems.append(f"{s}: {q['zero_frac']:.1%} of guides at zero "
                                f"(max {args.max_zero_frac:.0%}): dropout/"
                                f"bottleneck in a baseline arm")
            if not math.isnan(q["gini"]) and q["gini"] > args.max_gini:
                problems.append(f"{s}: Gini {q['gini']:.3f} (max {args.max_gini}): "
                                f"representation too skewed in a baseline arm")

    # replicate structure within control arms (swap detection)
    log = {s: np.log1p(mat[s]) for s in samples}
    by_cond = {}
    for s in samples:
        by_cond.setdefault(cond.get(s), []).append(s)

    rho_report = []
    for c in sorted(control_conditions):
        mates = by_cond.get(c, [])
        if len(mates) < 2:
            continue
        for s in mates:
            within = [spearman(log[s], log[o]) for o in mates if o != s]
            med_within = float(np.median(within))
            others = [spearman(log[s], log[o]) for o in samples
                      if cond.get(o) != c]
            best_other = max(others) if others else 0.0
            rho_report.append((s, c, med_within, best_other))
            if med_within < args.min_replicate_rho:
                problems.append(f"{s}: median Spearman to its {c} replicates is "
                                f"{med_within:.3f} (min {args.min_replicate_rho}): "
                                f"a mislabelled or swapped baseline library")
            if others and best_other > med_within:
                problems.append(f"{s}: correlates better with another condition "
                                f"({best_other:.3f}) than with its own {c} "
                                f"replicates ({med_within:.3f}): likely swapped")

    with open(args.out, "w") as fh:
        fh.write("sample\tcondition\treads\tmapped_frac\tcoverage\t"
                 "zero_frac\tgini\n")
        for s in samples:
            q = qc.get(s, {})
            fh.write(f"{s}\t{cond.get(s)}\t{q.get('reads', 0):.0f}\t"
                     f"{q.get('mapped_frac', 0):.4f}\t{q.get('coverage', 0):.1f}\t"
                     f"{q.get('zero_frac', 0):.4f}\t{q.get('gini', float('nan')):.4f}\n")
        if rho_report:
            fh.write("\n# control-arm replicate correlation\n")
            fh.write("# sample\tcondition\tmedian_rho_within\tbest_rho_other\n")
            for s, c, w, o in rho_report:
                fh.write(f"# {s}\t{c}\t{w:.4f}\t{o:.4f}\n")

    for s in samples:
        q = qc.get(s, {})
        print(f"  {s:<24} {str(cond.get(s)):<12} map {q.get('mapped_frac', 0):6.1%} "
              f"cov {q.get('coverage', 0):7.0f} zero {q.get('zero_frac', 0):5.1%} "
              f"gini {q.get('gini', float('nan')):.3f}")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
