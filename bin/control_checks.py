#!/usr/bin/env python3
"""
Did the screen behave? Two controls on a MAGeCK gene_summary, each a hard stop.

A gene_summary.txt is a ranked list with p-values and FDRs. It is produced
whether or not the assay worked and whether or not the null is calibrated, and
a mis-calibrated screen produces a longer, more exciting hit list, not an
error. Two internal controls say whether the numbers mean anything:

  negative control -- nontargeting guides.
    They target nothing, so they must behave like the null: roughly uniform
    p-values, and essentially none passing FDR. If controls look like hits, the
    normalisation or the counts are off and every FDR in the file is optimistic.

    Which file this is read from depends on how `mageck test` ran. Under
    --control-sgrna (the default here: control-normalised, control null) MAGeCK
    uses the control guides as the null and skips them as genes -- it logs
    "Skipping gene <id> for permutation" -- so they are absent from
    gene_summary.txt while each one is still scored in sgrna_summary.txt. That
    per-guide test is then the one used, and it is the better one: every
    control guide rather than a handful of pseudo-genes, from the same run
    being judged. If controls do appear as genes (a run without
    --control-sgrna), the gene-level test is used instead.

  positive control -- known essential genes (optional).
    In a proliferation or survival screen, guides against core-essential genes
    drop out over time. If a curated essential set does NOT deplete, the screen
    did not impose the selection you think it did, and a clean-looking enrichment
    list is being read off a screen that never worked.

Neither control looks at the genes of interest, so passing them does not make a
hit real; failing them means no hit can be.
"""
import argparse
import csv
import sys


def read_gene_summary(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def col(row, *names, default="nan"):
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return default


def load_list(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                out.add(s.split()[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-summary", required=True)
    ap.add_argument("--sgrna-summary", default=None,
                    help="sgrna_summary.txt from the same run; the control null "
                         "is read here when --control-sgrna kept the controls "
                         "out of gene_summary.txt")
    ap.add_argument("--control-sgrna-list", default=None,
                    help="file of exact control sgRNA ids, for --sgrna-summary")
    ap.add_argument("--control-genes", default="Control,NonTargeting,NTC",
                    help="comma-separated gene-label prefixes for controls")
    ap.add_argument("--control-genes-list", default=None,
                    help="file of exact control gene ids (overrides prefixes)")
    ap.add_argument("--essential-genes", default=None,
                    help="file of known essential gene ids (positive control)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--max-control-hit-frac", type=float, default=0.05,
                    help="fraction of control genes allowed below FDR")
    ap.add_argument("--min-control-median-p", type=float, default=0.3,
                    help="control p-values must not be systematically small")
    ap.add_argument("--essential-direction", choices=("neg", "pos"), default="neg")
    ap.add_argument("--min-essential-hit-frac", type=float, default=0.3,
                    help="fraction of essential genes that must pass FDR in "
                         "--essential-direction")
    ap.add_argument("--out", default="control_checks.tsv")
    args = ap.parse_args()

    rows = read_gene_summary(args.gene_summary)
    if not rows:
        raise SystemExit("gene_summary has no rows")

    if args.control_genes_list:
        ctrl_ids = load_list(args.control_genes_list)

        def is_control(gid):
            return gid in ctrl_ids
    else:
        toks = [t.strip().lower() for t in args.control_genes.split(",") if t.strip()]

        def is_control(gid):
            g = gid.lower()
            return any(g.startswith(t) for t in toks)

    problems = []
    lines = []

    # Controls as genes (a run without --control-sgrna), else per guide.
    control_rows = [r for r in rows if is_control(r["id"])]
    unit = "genes"
    ctrl_ps, ctrl_hits = [], []

    if control_rows:
        for r in control_rows:
            # Two-sided, so the floor below means the same thing here as it
            # does for sgrna_summary's p.twosided: min(neg,pos) alone is
            # uniform on [0, 0.5] under the null and would trip any floor
            # above 0.25 on a perfectly calibrated screen.
            one_sided = min(float(col(r, "neg|p-value", default="1")),
                            float(col(r, "pos|p-value", default="1")))
            ctrl_ps.append(min(1.0, 2.0 * one_sided))
            if (float(col(r, "neg|fdr")) < args.fdr or
                    float(col(r, "pos|fdr")) < args.fdr):
                ctrl_hits.append(r["id"])
    elif args.sgrna_summary:
        # MAGeCK skipped the controls as genes because they were the null;
        # every control guide is still scored in sgrna_summary.txt.
        unit = "sgRNAs"
        sg_rows = read_gene_summary(args.sgrna_summary)
        if args.control_sgrna_list:
            ctrl_sg = load_list(args.control_sgrna_list)

            def is_control_sg(row):
                return row.get("sgrna", "") in ctrl_sg
        else:
            def is_control_sg(row):
                return is_control(row.get("Gene", ""))

        for r in sg_rows:
            if not is_control_sg(r):
                continue
            ctrl_ps.append(float(col(r, "p.twosided", "p.low", default="1")))
            if float(col(r, "FDR", default="1")) < args.fdr:
                ctrl_hits.append(r.get("sgrna", "?"))

    n_ctrl = len(ctrl_ps)
    if n_ctrl == 0:
        problems.append("no controls found: they are not genes in gene_summary "
                        "and no --sgrna-summary gave per-guide scores, so the "
                        "null cannot be calibrated. Check --control_genes "
                        "matches the library's control labels.")
    else:
        frac = len(ctrl_hits) / n_ctrl
        lines.append(f"control {unit}: {n_ctrl}, below FDR {args.fdr}: "
                     f"{len(ctrl_hits)} ({frac:.1%})")
        if frac > args.max_control_hit_frac:
            problems.append(f"{frac:.1%} of control {unit} pass FDR {args.fdr} "
                            f"(max {args.max_control_hit_frac:.0%}), "
                            f"e.g. {ctrl_hits[:5]}: the null is not calibrated")

        # the control p-value distribution must not be skewed small
        ctrl_ps.sort()
        med = ctrl_ps[n_ctrl // 2]
        lines.append(f"control {unit} median p-value: {med:.3f}")
        if med < args.min_control_median_p:
            problems.append(f"control {unit} have median p-value {med:.3f} "
                            f"(min {args.min_control_median_p}): the control "
                            f"distribution is shifted, not uniform")

    if args.essential_genes:
        ess = load_list(args.essential_genes)
        d = args.essential_direction
        ess_rows = [r for r in rows if r["id"] in ess]
        found = len(ess_rows)
        if found == 0:
            problems.append(f"none of the {len(ess)} essential genes are in "
                            f"gene_summary; check the gene ids match the library")
        else:
            hits = [r["id"] for r in ess_rows if float(col(r, f"{d}|fdr")) < args.fdr]
            frac = len(hits) / found
            lines.append(f"essential genes: {found}/{len(ess)} present, "
                         f"{d}-selected below FDR: {len(hits)} ({frac:.1%})")
            if frac < args.min_essential_hit_frac:
                problems.append(f"only {frac:.1%} of essential genes are "
                                f"{d}-selected at FDR {args.fdr} "
                                f"(min {args.min_essential_hit_frac:.0%}): the "
                                f"selection did not take; the screen did not work")

    with open(args.out, "w") as fh:
        for ln in lines:
            fh.write(ln + "\n")

    for ln in lines:
        print("  " + ln)

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
