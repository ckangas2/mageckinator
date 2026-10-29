#!/usr/bin/env python3
"""
sgRNA library sanity, before a single read is counted.

MAGeCK aggregates sgRNA counts to a gene-level score, so the library file is
the map from guide to gene that the whole result rests on. A malformed library
does not error; it quietly mis-assigns reads and the gene ranking absorbs it:

  - two sgRNAs with the same spacer sequence but different gene labels: a read
    matching that sequence is counted once, for whichever row MAGeCK saw first.
    The two genes then share evidence, and both ranks are wrong.
  - a guide whose sequence also appears as a nontargeting control: real on-
    target signal leaks into the null, flattening it.
  - genes with far fewer guides than the rest: alpha-RRA ranks on the best few
    guides per gene, so a gene with one guide is compared on one data point and
    is neither reliably a hit nor reliably not.
  - no nontargeting controls at all: MAGeCK's permutation null and the
    control-based normalisation have nothing to stand on, and `control_checks`
    downstream cannot calibrate the FDR.

Writes a cleaned, MAGeCK-ready library (library.validated.csv: sgRNA,sequence,
gene) and the control-guide list used later for normalisation and the null.
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

DNA = re.compile(r"^[ACGT]+$")


def read_library(path):
    with open(path, newline="") as fh:
        sniff = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sniff and "," not in sniff else ","
        rows = list(csv.reader(fh, delimiter=delim))
    if not rows:
        raise SystemExit(f"{path} is empty")
    # tolerate a header line (first row not DNA in column 2)
    start = 0
    if len(rows[0]) >= 2 and not DNA.match(rows[0][1].strip().upper()):
        start = 1
    out = []
    for i, r in enumerate(rows[start:], start=start + 1):
        if not r or all(not c.strip() for c in r):
            continue
        if len(r) < 3:
            raise SystemExit(f"{path} line {i}: expected sgRNA,sequence,gene; "
                             f"got {r}")
        out.append((r[0].strip(), r[1].strip().upper(), r[2].strip()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True)
    ap.add_argument("--control-genes", default="Control,NonTargeting,NTC",
                    help="comma-separated gene labels that mark nontargeting "
                         "controls (case-insensitive, prefix match)")
    ap.add_argument("--min-guides", type=int, default=3,
                    help="minimum sgRNAs per targeted gene")
    ap.add_argument("--require-controls", action="store_true", default=True)
    ap.add_argument("--allow-no-controls", dest="require_controls",
                    action="store_false")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    rows = read_library(args.library)
    problems = []

    ctrl_tokens = [t.strip().lower() for t in args.control_genes.split(",") if t.strip()]

    def is_control(gene):
        g = gene.lower()
        return any(g.startswith(t) for t in ctrl_tokens)

    # unique sgRNA ids
    id_counts = Counter(r[0] for r in rows)
    dup_ids = [k for k, v in id_counts.items() if v > 1]
    if dup_ids:
        problems.append(f"{len(dup_ids)} sgRNA id(s) appear more than once, "
                        f"e.g. {dup_ids[:5]}")

    # sequence validity and uniform length
    lengths = Counter(len(r[1]) for r in rows)
    bad_seq = [r[0] for r in rows if not DNA.match(r[1])]
    if bad_seq:
        problems.append(f"{len(bad_seq)} sgRNA(s) have non-ACGT sequence, "
                        f"e.g. {bad_seq[:5]}")
    if len(lengths) > 1:
        problems.append(f"sgRNA sequences are not a uniform length: "
                        f"{dict(lengths)} (mageck count needs one --sgrna-len)")

    # one sequence -> one gene
    seq_to_genes = {}
    for sg, seq, gene in rows:
        seq_to_genes.setdefault(seq, set()).add(gene)
    collisions = {seq: g for seq, g in seq_to_genes.items() if len(g) > 1}
    if collisions:
        ex = list(collisions.items())[:5]
        problems.append(f"{len(collisions)} spacer sequence(s) are shared across "
                        f"genes; a read cannot be assigned unambiguously, "
                        f"e.g. " + "; ".join(f"{s}->{sorted(g)}" for s, g in ex))

    # a control sequence that is also a targeting guide
    control_seqs = {seq for sg, seq, gene in rows if is_control(gene)}
    target_seqs = {seq for sg, seq, gene in rows if not is_control(gene)}
    leaked = control_seqs & target_seqs
    if leaked:
        problems.append(f"{len(leaked)} control spacer(s) are also used as "
                        f"targeting guides; the null is contaminated with "
                        f"on-target signal")

    # guides per targeted gene
    guides_per_gene = Counter(gene for sg, seq, gene in rows if not is_control(gene))
    thin = {g: n for g, n in guides_per_gene.items() if n < args.min_guides}
    if thin:
        ex = sorted(thin.items())[:5]
        problems.append(f"{len(thin)} gene(s) have < {args.min_guides} guides, "
                        f"e.g. {ex}")

    # controls present
    n_controls = sum(1 for sg, seq, gene in rows if is_control(gene))
    if n_controls == 0 and args.require_controls:
        problems.append("no nontargeting control guides found "
                        f"(looked for genes starting with {ctrl_tokens}); "
                        "pass --control-genes or --allow-no-controls")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1

    with open(out_dir / "library.validated.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sgRNA", "sequence", "gene"])
        w.writerows(rows)

    with open(out_dir / "control_sgrnas.txt", "w") as fh:
        for sg, seq, gene in rows:
            if is_control(gene):
                fh.write(sg + "\n")

    with open(out_dir / "control_genes.txt", "w") as fh:
        for gene in sorted({gene for sg, seq, gene in rows if is_control(gene)}):
            fh.write(gene + "\n")

    n_genes = len(guides_per_gene)
    seq_len = next(iter(lengths))
    with open(out_dir / "library_report.txt", "w") as fh:
        fh.write("mageckinator library report\n")
        fh.write(f"library        : {args.library}\n")
        fh.write(f"sgRNAs         : {len(rows)}\n")
        fh.write(f"targeted genes : {n_genes}\n")
        fh.write(f"control sgRNAs : {n_controls}\n")
        fh.write(f"spacer length  : {seq_len}\n")
        med = sorted(guides_per_gene.values())[n_genes // 2] if n_genes else 0
        fh.write(f"median guides/gene : {med}\n")

    print(f"  {len(rows)} sgRNAs, {n_genes} targeted genes, "
          f"{n_controls} controls, spacer {seq_len} bp  ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
