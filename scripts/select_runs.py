#!/usr/bin/env python3
"""
Turn an ENA run report into a samplesheet + download list for the GSE158298
pooled CRISPR screen, keeping only the screen libraries and labelling their
condition from the sample title.

    python3 scripts/select_runs.py ena_runs.tsv --list
    python3 scripts/select_runs.py ena_runs.tsv --fastq-dir fastq

Condition is parsed, not hand-copied: a title mentioning "plasmid" or
"library" is the Plasmid baseline; "low" + "MOI" is LowMOI; "high" + "MOI" is
HighMOI. Anything that does not match (the study's scRNA-seq, ECCITE-seq and
bulk RNA-seq runs) is dropped. Inspect --list first; if the titles do not
parse as expected, fix the patterns here rather than editing the samplesheet
by hand, so the selection stays reproducible.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict


def classify(title):
    t = title.lower()
    is_moi = "moi" in t or "infect" in t
    if "plasmid" in t or "library" in t or ("pdna" in t):
        return "Plasmid"
    if is_moi and ("low" in t or "0.01" in t):
        return "LowMOI"
    if is_moi and ("high" in t or "0.3" in t):
        return "HighMOI"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report")
    ap.add_argument("--fastq-dir", default="fastq")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    with open(args.report, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    kept = []
    counters = defaultdict(int)
    for r in rows:
        title = r.get("sample_title") or r.get("sample_alias") or r.get("run_accession")
        strat = (r.get("library_strategy") or "").upper()
        cond = classify(title)
        # a pooled CRISPR screen reads as amplicon/OTHER, never RNA-Seq
        if cond and strat not in ("RNA-SEQ", "SINGLE CELL"):
            counters[cond] += 1
            kept.append((r, cond, counters[cond]))

    if args.list:
        print("run_accession\tcondition\treplicate\tstrategy\tsample_title")
        for r, cond, rep in kept:
            print(f"{r['run_accession']}\t{cond}\t{rep}\t"
                  f"{r.get('library_strategy','')}\t{r.get('sample_title','')}")
        # also show what was dropped, so a mis-parse is visible
        for r in rows:
            if not any(r is k[0] for k in kept):
                print(f"{r['run_accession']}\t-\t-\t"
                      f"{r.get('library_strategy','')}\t{r.get('sample_title','')}")
        return 0

    if not kept:
        sys.exit("no screen libraries matched; run with --list and adjust classify()")

    downloads = []
    with open("samplesheet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "fastq_1", "condition", "replicate"])
        for r, cond, rep in kept:
            sample = f"{cond}_rep{rep}"
            # screen libraries are single-end amplicon; take the first fastq
            urls = (r.get("fastq_ftp") or "").split(";")
            md5s = (r.get("fastq_md5") or "").split(";")
            if not urls or not urls[0]:
                sys.exit(f"{r['run_accession']}: no fastq_ftp in the ENA report")
            url = urls[0] if urls[0].startswith("http") else "https://" + urls[0]
            fq = f"{args.fastq_dir}/{r['run_accession']}.fastq.gz"
            w.writerow([sample, fq, cond, rep])
            downloads.append((url, md5s[0], fq))

    with open("config/downloads.tsv", "w", newline="") as fh:
        for url, md5, fq in downloads:
            fh.write(f"{url}\t{md5}\t{fq}\n")

    for cond in ("Plasmid", "LowMOI", "HighMOI"):
        print(f"  {cond}: {counters[cond]} library(ies)")
    print(f"selected {len(kept)} runs -> samplesheet.csv, config/downloads.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
