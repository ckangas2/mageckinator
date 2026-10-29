#!/usr/bin/env python3
"""
Unit tests: one planted failure per gate, plus the passing case.

    bash tests/run_tests.sh        (or: python3 -m pytest tests -q)

No MAGeCK, no Docker, no network. Each bin/ script is run as a subprocess,
exactly as the pipeline runs it, and judged on its exit code (0 = gate passed,
non-zero = gate stopped the run) and the files it writes. The MAGeCK-format
tables the gates consume are built inline here so the tests need no tools; the
end-to-end recovery of the planted truth is checked separately, against a real
MAGeCK run, by check_test_results.py.
"""
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
DATA = ROOT / "tests" / "data"


def run(script, *args, cwd=None):
    p = subprocess.run([sys.executable, "-I", str(BIN / script), *map(str, args)],
                       capture_output=True, text=True, cwd=cwd)
    return p.returncode, p.stdout + p.stderr


def write(path, text):
    Path(path).write_text(text)
    return path


# ------------------------------------------------------------- fixtures

def samplesheet(tmp, rows, name="ss.csv"):
    hdr = "sample,fastq_1,condition,replicate\n"
    body = "".join(f"{s},{s}.fq.gz,{c},{r}\n" for s, c, r in rows)
    return write(tmp / name, hdr + body)


def comparisons(tmp, body, name="cmp.csv"):
    return write(tmp / name, "id,control,treatment,selection\n" + body)


def library(tmp, rows, name="lib.csv"):
    body = "".join(f"{sg},{seq},{gene}\n" for sg, seq, gene in rows)
    return write(tmp / name, "sgRNA,sequence,gene\n" + body)


GS_COLS = ["id", "num",
           "neg|score", "neg|p-value", "neg|fdr", "neg|rank", "neg|goodsgrna", "neg|lfc",
           "pos|score", "pos|p-value", "pos|fdr", "pos|rank", "pos|goodsgrna", "pos|lfc"]


def gene_summary(tmp, genes, name="gs.txt"):
    """genes: dict id -> dict of overrides for the columns above."""
    lines = ["\t".join(GS_COLS)]
    for i, (gid, ov) in enumerate(genes.items(), 1):
        row = {c: "1" for c in GS_COLS}
        row.update({"id": gid, "num": "4",
                    "neg|rank": str(i), "pos|rank": str(i),
                    "neg|lfc": "0", "pos|lfc": "0",
                    "neg|goodsgrna": "0", "pos|goodsgrna": "0"})
        row.update({k: str(v) for k, v in ov.items()})
        lines.append("\t".join(row[c] for c in GS_COLS))
    return write(tmp / name, "\n".join(lines) + "\n")


def count_tables(tmp, samples, per_gene):
    """Build counts.txt + countsummary.txt from a {sample: {sgRNA: count}}."""
    sgrnas = sorted(next(iter(per_gene.values())).keys())
    gene_of = {sg: sg.split("_sg")[0] for sg in sgrnas}
    # counts.txt
    lines = ["sgRNA\tgene\t" + "\t".join(samples)]
    for sg in sgrnas:
        lines.append(sg + "\t" + gene_of[sg] + "\t" +
                     "\t".join(str(per_gene[s][sg]) for s in samples))
    counts = write(tmp / "counts.txt", "\n".join(lines) + "\n")
    # countsummary.txt
    cs = ["\t".join(["File", "Label", "Reads", "Mapped", "Percentage",
                     "TotalsgRNAs", "Zerocounts", "GiniIndex"])]
    for s in samples:
        mapped = sum(per_gene[s].values())
        zero = sum(1 for v in per_gene[s].values() if v == 0)
        cs.append("\t".join([f"{s}.fq.gz", s, str(int(mapped / 0.95)), str(mapped),
                             "0.95", str(len(sgrnas)), str(zero), "0.10"]))
    summary = write(tmp / "countsummary.txt", "\n".join(cs) + "\n")
    return counts, summary


# =================================================================== check_design

def test_design_passes_on_the_test_dataset(tmp_path):
    rc, out = run("check_design.py",
                  "--samplesheet", DATA / "samplesheet.csv",
                  "--comparisons", DATA / "comparisons.csv",
                  "--out-dir", tmp_path)
    assert rc == 0, out
    assert (tmp_path / "comparisons.validated.tsv").exists()


def test_design_fails_when_control_equals_treatment(tmp_path):
    ss = samplesheet(tmp_path, [("P1", "Plasmid", "1"), ("P2", "Plasmid", "2"),
                                ("H1", "HighMOI", "1"), ("H2", "HighMOI", "2")])
    cmp = comparisons(tmp_path, "c1,Plasmid,Plasmid,positive\n")
    rc, out = run("check_design.py", "--samplesheet", ss, "--comparisons", cmp,
                  "--out-dir", tmp_path)
    assert rc != 0 and "two different conditions" in out


def test_design_fails_on_unknown_condition(tmp_path):
    ss = samplesheet(tmp_path, [("P1", "Plasmid", "1"), ("P2", "Plasmid", "2"),
                                ("H1", "HighMOI", "1"), ("H2", "HighMOI", "2")])
    cmp = comparisons(tmp_path, "c1,Plasmid,LowMOI,positive\n")
    rc, out = run("check_design.py", "--samplesheet", ss, "--comparisons", cmp,
                  "--out-dir", tmp_path)
    assert rc != 0 and "no samples" in out


def test_design_fails_when_an_arm_is_underpowered(tmp_path):
    ss = samplesheet(tmp_path, [("P1", "Plasmid", "1"), ("P2", "Plasmid", "2"),
                                ("H1", "HighMOI", "1")])
    cmp = comparisons(tmp_path, "c1,Plasmid,HighMOI,positive\n")
    rc, out = run("check_design.py", "--samplesheet", ss, "--comparisons", cmp,
                  "--min-replicates", "2", "--out-dir", tmp_path)
    assert rc != 0 and "need 2" in out


# =================================================================== validate_library

def test_library_passes_on_the_test_dataset(tmp_path):
    rc, out = run("validate_library.py", "--library", DATA / "library.csv",
                  "--out-dir", tmp_path)
    assert rc == 0, out
    assert (tmp_path / "library.validated.csv").exists()
    assert (tmp_path / "control_sgrnas.txt").read_text().strip()


def test_library_fails_on_sequence_shared_across_genes(tmp_path):
    lib = library(tmp_path, [
        ("GENE1_sg1", "ACGTACGTACGTACGTACGT", "GENE1"),
        ("GENE1_sg2", "TTTTACGTACGTACGTACGT", "GENE1"),
        ("GENE1_sg3", "CCCCACGTACGTACGTACGT", "GENE1"),
        ("GENE2_sg1", "ACGTACGTACGTACGTACGT", "GENE2"),  # same spacer, other gene
        ("GENE2_sg2", "GGGGACGTACGTACGTACGT", "GENE2"),
        ("GENE2_sg3", "AAAAACGTACGTACGTACGT", "GENE2"),
        ("Control_1_sg1", "GAGAGAGAGAGAGAGAGAGA", "Control_1"),
    ])
    rc, out = run("validate_library.py", "--library", lib, "--min-guides", "3",
                  "--out-dir", tmp_path)
    assert rc != 0 and "shared across" in out


def test_library_fails_on_too_few_guides_per_gene(tmp_path):
    lib = library(tmp_path, [
        ("GENE1_sg1", "ACGTACGTACGTACGTAC01".replace("0", "A").replace("1", "C"), "GENE1"),
        ("GENE1_sg2", "TTTTACGTACGTACGTACGT", "GENE1"),
        ("Control_1_sg1", "GAGAGAGAGAGAGAGAGAGA", "Control_1"),
    ])
    rc, out = run("validate_library.py", "--library", lib, "--min-guides", "3",
                  "--out-dir", tmp_path)
    assert rc != 0 and "guides" in out


def test_library_fails_when_no_controls(tmp_path):
    lib = library(tmp_path, [
        ("GENE1_sg1", "ACGTACGTACGTACGTACGA", "GENE1"),
        ("GENE1_sg2", "TTTTACGTACGTACGTACGT", "GENE1"),
        ("GENE1_sg3", "CCCCACGTACGTACGTACGT", "GENE1"),
    ])
    rc, out = run("validate_library.py", "--library", lib, "--min-guides", "3",
                  "--out-dir", tmp_path)
    assert rc != 0 and "control" in out.lower()


# =================================================================== count_qc

def _good_counts(samples=("Plasmid_rep1", "Plasmid_rep2",
                          "HighMOI_rep1", "HighMOI_rep2", "HighMOI_rep3")):
    # 20 sgRNAs (5 genes x 4); host factor GENE1 enriched in HighMOI
    per = {s: {} for s in samples}
    base = {}
    for g in range(1, 6):
        for k in range(1, 5):
            base[f"GENE{g}_sg{k}"] = 50 + g * 3 + k
    for s in samples:
        for sg, b in base.items():
            v = b
            if s.startswith("HighMOI") and sg.startswith("GENE1_"):
                v = b * 5
            per[s][sg] = v
    return samples, per


def test_count_qc_passes_on_clean_counts(tmp_path):
    samples, per = _good_counts()
    counts, summ = count_tables(tmp_path, list(samples), per)
    ss = samplesheet(tmp_path, [(s, s.split("_rep")[0], s[-1]) for s in samples])
    rc0, o0 = run("check_design.py", "--samplesheet", ss,
                  "--comparisons", comparisons(tmp_path, "c1,Plasmid,HighMOI,positive\n"),
                  "--out-dir", tmp_path)
    assert rc0 == 0, o0
    rc, out = run("count_qc.py", "--countsummary", summ, "--counts", counts,
                  "--samplesheet", ss,
                  "--comparisons", tmp_path / "comparisons.validated.tsv",
                  "--min-coverage", "50", "--out", tmp_path / "qc.tsv")
    assert rc == 0, out


def test_count_qc_fails_on_low_mapping(tmp_path):
    samples, per = _good_counts()
    counts, summ = count_tables(tmp_path, list(samples), per)
    # clobber one sample's mapping rate in the countsummary
    text = summ.read_text().splitlines()
    text[1] = text[1].replace("\t0.95\t", "\t0.20\t")
    summ.write_text("\n".join(text) + "\n")
    ss = samplesheet(tmp_path, [(s, s.split("_rep")[0], s[-1]) for s in samples])
    run("check_design.py", "--samplesheet", ss,
        "--comparisons", comparisons(tmp_path, "c1,Plasmid,HighMOI,positive\n"),
        "--out-dir", tmp_path)
    rc, out = run("count_qc.py", "--countsummary", summ, "--counts", counts,
                  "--samplesheet", ss,
                  "--comparisons", tmp_path / "comparisons.validated.tsv",
                  "--min-coverage", "50", "--out", tmp_path / "qc.tsv")
    assert rc != 0 and "map to the library" in out


def test_count_qc_fails_on_swapped_baseline(tmp_path):
    samples, per = _good_counts()
    # Make Plasmid_rep2 actually look like a HighMOI survivor pool: swap its
    # profile. Its correlation to the other Plasmid drops below the floor.
    per["Plasmid_rep2"] = dict(per["HighMOI_rep1"])
    counts, summ = count_tables(tmp_path, list(samples), per)
    ss = samplesheet(tmp_path, [(s, s.split("_rep")[0], s[-1]) for s in samples])
    run("check_design.py", "--samplesheet", ss,
        "--comparisons", comparisons(tmp_path, "c1,Plasmid,HighMOI,positive\n"),
        "--out-dir", tmp_path)
    rc, out = run("count_qc.py", "--countsummary", summ, "--counts", counts,
                  "--samplesheet", ss,
                  "--comparisons", tmp_path / "comparisons.validated.tsv",
                  "--min-coverage", "50", "--min-replicate-rho", "0.8",
                  "--out", tmp_path / "qc.tsv")
    assert rc != 0 and ("swapped" in out or "replicates" in out)


# =================================================================== control_checks

def _calibrated_gene_summary(tmp_path):
    genes = {}
    # host factors: strong positive selection
    for g in ["ACE2", "RAB7A", "CTSL"]:
        genes[g] = {"pos|p-value": "1e-6", "pos|fdr": "0.001", "pos|lfc": "2.3"}
    # essentials: strong negative selection
    for g in ["RPL9", "RPS19", "POLR2A"]:
        genes[g] = {"neg|p-value": "1e-6", "neg|fdr": "0.001", "neg|lfc": "-3.0"}
    # flat genes: null
    for i in range(1, 6):
        genes[f"GENE{i:03d}"] = {"neg|p-value": "0.5", "pos|p-value": "0.5"}
    # control pseudo-genes: null, uniform-ish p-values
    for i, p in enumerate([0.3, 0.45, 0.6, 0.72, 0.9, 0.25, 0.55, 0.8], 1):
        genes[f"Control_{i:02d}"] = {"neg|p-value": str(p), "pos|p-value": str(p),
                                     "neg|fdr": "0.9", "pos|fdr": "0.9"}
    return gene_summary(tmp_path, genes)


def test_control_checks_pass_when_calibrated(tmp_path):
    gs = _calibrated_gene_summary(tmp_path)
    ess = write(tmp_path / "ess.txt", "RPL9\nRPS19\nPOLR2A\n")
    rc, out = run("control_checks.py", "--gene-summary", gs,
                  "--essential-genes", ess, "--out", tmp_path / "cc.tsv")
    assert rc == 0, out


def test_control_checks_fail_when_controls_are_hits(tmp_path):
    genes = {"ACE2": {"pos|fdr": "0.001"}}
    # half the control pseudo-genes pass FDR: the null is hot
    for i in range(1, 9):
        fdr = "0.001" if i <= 4 else "0.9"
        p = "1e-5" if i <= 4 else "0.6"
        genes[f"Control_{i:02d}"] = {"pos|fdr": fdr, "neg|fdr": "0.9",
                                     "pos|p-value": p, "neg|p-value": "0.6"}
    gs = gene_summary(tmp_path, genes)
    rc, out = run("control_checks.py", "--gene-summary", gs,
                  "--out", tmp_path / "cc.tsv")
    assert rc != 0 and ("not calibrated" in out or "shifted" in out)


def test_control_checks_fail_when_essentials_do_not_deplete(tmp_path):
    genes = {}
    for g in ["RPL9", "RPS19", "POLR2A"]:
        genes[g] = {"neg|fdr": "0.9", "neg|p-value": "0.7"}   # did NOT deplete
    for i in range(1, 9):
        genes[f"Control_{i:02d}"] = {"neg|p-value": str(0.3 + i * 0.05),
                                     "pos|p-value": str(0.3 + i * 0.05),
                                     "neg|fdr": "0.9", "pos|fdr": "0.9"}
    gs = gene_summary(tmp_path, genes)
    ess = write(tmp_path / "ess.txt", "RPL9\nRPS19\nPOLR2A\n")
    rc, out = run("control_checks.py", "--gene-summary", gs,
                  "--essential-genes", ess, "--out", tmp_path / "cc.tsv")
    assert rc != 0 and "did not take" in out


# =================================================================== reconcile

def test_reconcile_passes_when_hits_match_published(tmp_path):
    genes = {}
    for i, g in enumerate(["ACE2", "RAB7A", "CTSL", "VPS35"], 1):
        genes[g] = {"pos|rank": str(i), "pos|fdr": "0.001"}
    for i in range(1, 6):
        genes[f"GENE{i:03d}"] = {"pos|rank": str(10 + i), "pos|fdr": "0.9"}
    gs = gene_summary(tmp_path, genes)
    pub = write(tmp_path / "pub.tsv", "gene\trank\nACE2\t1\nRAB7A\t2\nCTSL\t3\nVPS35\t4\n")
    rc, out = run("reconcile.py", "--gene-summary", gs, "--published", pub,
                  "--selection", "positive", "--top-n", "5",
                  "--min-recall", "0.5", "--out", tmp_path / "rec.tsv")
    assert rc == 0, out


def test_reconcile_fails_when_direction_is_inverted(tmp_path):
    # published hits are actually enriched, but our file has them as the top
    # of the NEGATIVE list (as if treatment/control were swapped)
    genes = {}
    for i, g in enumerate(["ACE2", "RAB7A", "CTSL", "VPS35"], 1):
        genes[g] = {"neg|rank": str(i), "neg|fdr": "0.001",
                    "pos|rank": str(500 + i), "pos|fdr": "0.9"}
    for i in range(1, 6):
        genes[f"GENE{i:03d}"] = {"pos|rank": str(i), "pos|fdr": "0.5"}
    gs = gene_summary(tmp_path, genes)
    pub = write(tmp_path / "pub.tsv", "gene\nACE2\nRAB7A\nCTSL\nVPS35\n")
    rc, out = run("reconcile.py", "--gene-summary", gs, "--published", pub,
                  "--selection", "positive", "--top-n", "5",
                  "--min-recall", "0.5", "--out", tmp_path / "rec.tsv")
    assert rc != 0 and ("swapped" in out or "not the same screen" in out)


# =================================================================== collect_summary

def test_collect_summary_runs(tmp_path):
    samples, per = _good_counts()
    counts, summ = count_tables(tmp_path, list(samples), per)
    ss = samplesheet(tmp_path, [(s, s.split("_rep")[0], s[-1]) for s in samples])
    run("check_design.py", "--samplesheet", ss,
        "--comparisons", comparisons(tmp_path, "c1,Plasmid,HighMOI,positive\n"),
        "--out-dir", tmp_path)
    run("count_qc.py", "--countsummary", summ, "--counts", counts,
        "--samplesheet", ss, "--comparisons", tmp_path / "comparisons.validated.tsv",
        "--min-coverage", "50", "--out", tmp_path / "count_qc.tsv")
    gs = _calibrated_gene_summary(tmp_path)
    rc, out = run("collect_summary.py", "--count-qc", tmp_path / "count_qc.tsv",
                  "--comparisons", tmp_path / "comparisons.validated.tsv",
                  "--gene-summary", f"c1={gs}", "--out-dir", tmp_path / "sum")
    assert rc == 0, out
    assert (tmp_path / "sum" / "comparisons.tsv").exists()
