#!/usr/bin/env python3
"""
Derive a MAGeCK-style countsummary.txt from an existing count table.

When a run starts from a deposited count matrix (sgRNA, gene, <samples>)
instead of FASTQ, there is no `mageck count` summary. count_qc still wants
one, so this computes the per-sample fields that can be known from counts
alone: total reads, zero-count guides, Gini index, coverage. Mapping rate is
unknowable from a table, so it is reported as 1.0 and the mapping gate is a
no-op on this path -- stated plainly rather than faked to a plausible 95%.
"""
import argparse
import sys

import numpy as np


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    cum = np.cumsum(x)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out", default="countsummary.txt")
    args = ap.parse_args()

    with open(args.counts) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[2:]
        cols = {s: [] for s in samples}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for s, v in zip(samples, f[2:]):
                cols[s].append(float(v))

    with open(args.out, "w") as fh:
        fh.write("File\tLabel\tReads\tMapped\tPercentage\tTotalsgRNAs\t"
                 "Zerocounts\tGiniIndex\n")
        for s in samples:
            v = np.asarray(cols[s])
            mapped = int(v.sum())
            fh.write(f"{s}.table\t{s}\t{mapped}\t{mapped}\t1.0\t{len(v)}\t"
                     f"{int((v == 0).sum())}\t{gini(v):.4f}\n")
    print(f"  wrote countsummary for {len(samples)} samples from {args.counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
