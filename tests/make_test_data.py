#!/usr/bin/env python3
"""
Regenerate tests/data/ for the test profile, deterministically.

    python3 tests/make_test_data.py            # writes tests/data/
    python3 tests/make_test_data.py regen      # writes ./regen/ (CI compares)

A small pooled CRISPR survival screen, planted so the answer is known:

  library : 50 targeted genes x 4 guides + 20 nontargeting control
            pseudo-genes x 4 guides  (= 280 sgRNAs, 20 bp spacers)
  design  : condition Plasmid (2 libraries), condition HighMOI (3 libraries)
  biology : HighMOI is the survivor pool of a killing screen.
            - 8 host-factor genes are ENRICHED ~5x in survivors (positive
              selection): knocking them out blocks infection, so those cells
              live. These are the hits, and they carry real SARS-CoV-2 screen
              names so published_hits.tsv is meaningful.
            - 8 essential genes DEPLETE ~8x (negative selection): the assay's
              positive control; if they do not drop, the screen did not work.
            - 34 genes and all controls stay flat.

check_test_results.py asserts, against real `mageck count` + `mageck test`
output, that the host factors come back as the top positive-selection hits,
the essentials as depleted, and the controls as null.

numpy's PCG64 is stable across versions, so the reads regenerate byte-for-byte
(CI checks this on decompressed content).
"""
import gzip
import sys
from pathlib import Path

import numpy as np

SEED = 20201024  # GSE158298 deposit era
RNG = np.random.default_rng(SEED)

# 5' context ends in the canonical ...CACCG; 3' is the tracrRNA start. A real
# vector read is longer, but this is the part `mageck count` keys on.
ADAPTER5 = "CACCG"
ADAPTER3 = "GTTTAGAGC"
SPACER_LEN = 20

HOST_FACTORS = ["ACE2", "RAB7A", "CTSL", "ATP6V1A", "ATP6AP1",
                "VPS35", "NPC1", "PIK3C3"]
ESSENTIALS = ["RPL9", "RPS19", "POLR2A", "PSMB3", "EIF3B",
              "CDK1", "RRM1", "SNRPD1"]
N_FLAT = 34
N_CONTROL_GENES = 20
GUIDES_PER_GENE = 4

PLASMID_REPS = 2
HIGHMOI_REPS = 3
BASE_LOW, BASE_HIGH = 30, 70     # per-guide library representation
ENRICH = 5.0                     # host factors in survivors
DEPLETE = 8.0                    # essentials in survivors


def rand_spacer(rng):
    return "".join(rng.choice(list("ACGT"), size=SPACER_LEN))


def build_library(rng):
    """Return list of (sgRNA_id, spacer, gene, role) and a used-sequence set."""
    rows = []
    used = set()

    def guides(gene, role, n=GUIDES_PER_GENE):
        for i in range(1, n + 1):
            while True:
                s = rand_spacer(rng)
                if s not in used:
                    used.add(s)
                    break
            rows.append((f"{gene}_sg{i}", s, gene, role))

    for g in HOST_FACTORS:
        guides(g, "host_factor")
    for g in ESSENTIALS:
        guides(g, "essential")
    for i in range(1, N_FLAT + 1):
        guides(f"GENE{i:03d}", "flat")
    for i in range(1, N_CONTROL_GENES + 1):
        guides(f"Control_{i:02d}", "control")
    return rows


def planted_counts(rows, rng):
    """Per-guide integer counts for every sample, as a dict sample->list."""
    base = {sg: int(rng.integers(BASE_LOW, BASE_HIGH + 1)) for sg, _, _, _ in rows}

    samples = ([f"Plasmid_rep{i}" for i in range(1, PLASMID_REPS + 1)] +
               [f"HighMOI_rep{i}" for i in range(1, HIGHMOI_REPS + 1)])
    counts = {s: [] for s in samples}
    for sg, _spacer, _gene, role in rows:
        b = base[sg]
        for s in samples:
            if s.startswith("Plasmid"):
                lam = b
            else:
                if role == "host_factor":
                    lam = b * ENRICH
                elif role == "essential":
                    lam = b / DEPLETE
                else:
                    lam = b
            counts[s].append(int(rng.poisson(max(lam, 0.5))))
    return samples, counts


def write_fastq(path, rows, sample_counts):
    """One read per counted molecule: ADAPTER5 + spacer + ADAPTER3."""
    lines = []
    read = 0
    for (sg, spacer, _gene, _role), n in zip(rows, sample_counts):
        seq = ADAPTER5 + spacer + ADAPTER3
        qual = "I" * len(seq)
        for _ in range(n):
            read += 1
            lines.append(f"@{path.stem}_{read}\n{seq}\n+\n{qual}\n")
    data = "".join(lines).encode()
    with gzip.GzipFile(path, "wb", mtime=0) as fh:   # mtime=0 -> reproducible
        fh.write(data)


def main():
    out = Path("regen" if len(sys.argv) > 1 and sys.argv[1] == "regen"
               else Path(__file__).resolve().parent / "data")
    out.mkdir(parents=True, exist_ok=True)

    rows = build_library(RNG)
    samples, counts = planted_counts(rows, RNG)

    # library.csv
    with open(out / "library.csv", "w") as fh:
        fh.write("sgRNA,sequence,gene\n")
        for sg, spacer, gene, _role in rows:
            fh.write(f"{sg},{spacer},{gene}\n")

    # control sgRNA list (for mageck --control-sgrna normalisation)
    with open(out / "control_sgrnas.txt", "w") as fh:
        for sg, _s, gene, role in rows:
            if role == "control":
                fh.write(sg + "\n")

    # FASTQ per sample
    for i, s in enumerate(samples):
        write_fastq(out / f"{s}.fastq.gz", rows, [counts[s][j] for j in range(len(rows))])

    # samplesheet.csv
    with open(out / "samplesheet.csv", "w") as fh:
        fh.write("sample,fastq_1,condition,replicate\n")
        for s in samples:
            cond, rep = s.split("_rep")
            fh.write(f"{s},{s}.fastq.gz,{cond},{rep}\n")

    # comparisons.csv
    with open(out / "comparisons.csv", "w") as fh:
        fh.write("id,control,treatment,selection\n")
        fh.write("hostfactors_highmoi_vs_plasmid,Plasmid,HighMOI,positive\n")

    # truth.tsv
    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene\trole\n")
        for gene, role in dict((g, r) for _sg, _s, g, r in rows).items():
            fh.write(f"{gene}\t{role}\n")

    # published_hits.tsv (what the "paper" reported: the host factors)
    with open(out / "published_hits.tsv", "w") as fh:
        fh.write("gene\trank\n")
        for i, g in enumerate(HOST_FACTORS, 1):
            fh.write(f"{g}\t{i}\n")

    # essential gene list (positive control)
    with open(out / "essential_genes.txt", "w") as fh:
        for g in ESSENTIALS:
            fh.write(g + "\n")

    # deposited count table (the "start from counts" fast path)
    with open(out / "deposited_counts.tsv", "w") as fh:
        fh.write("sgRNA\tgene\t" + "\t".join(samples) + "\n")
        for j, (sg, _spacer, gene, _role) in enumerate(rows):
            vals = "\t".join(str(counts[s][j]) for s in samples)
            fh.write(f"{sg}\t{gene}\t{vals}\n")

    # ---- negative-case fixtures: one planted failure per gate (CI) ----
    # check_design: treatment == control
    with open(out / "comparisons_bad.csv", "w") as fh:
        fh.write("id,control,treatment,selection\n")
        fh.write("self_vs_self,Plasmid,Plasmid,positive\n")

    # validate_library: a spacer shared across two genes
    lib_rows = [(sg, spacer, gene) for sg, spacer, gene, _r in rows]
    dup_spacer = lib_rows[0][1]            # reuse ACE2_sg1's spacer for a flat gene
    with open(out / "library_dup.csv", "w") as fh:
        fh.write("sgRNA,sequence,gene\n")
        for sg, spacer, gene in lib_rows:
            fh.write(f"{sg},{spacer},{gene}\n")
        fh.write(f"GENE001_sgX,{dup_spacer},GENE001\n")   # collides with ACE2_sg1

    # count_qc: a swapped baseline (Plasmid_rep2 column is really a survivor pool)
    with open(out / "counts_swapped.tsv", "w") as fh:
        fh.write("sgRNA\tgene\t" + "\t".join(samples) + "\n")
        for j, (sg, _spacer, gene, _role) in enumerate(rows):
            vals = []
            for s in samples:
                src = "HighMOI_rep1" if s == "Plasmid_rep2" else s
                vals.append(str(counts[src][j]))
            fh.write(f"{sg}\t{gene}\t" + "\t".join(vals) + "\n")

    # control_checks: an essential list that is actually flat genes -> no dropout
    with open(out / "fake_essentials.txt", "w") as fh:
        for i in range(1, 11):
            fh.write(f"GENE{i:03d}\n")

    # reconcile: a published list of flat genes -> recall collapses
    with open(out / "wrong_published.tsv", "w") as fh:
        fh.write("gene\n")
        for i in range(11, 21):
            fh.write(f"GENE{i:03d}\n")

    print(f"wrote {len(rows)} sgRNAs, {len(samples)} libraries to {out}/")


if __name__ == "__main__":
    main()
