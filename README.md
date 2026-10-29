# mageckinator

[![CI](https://github.com/ckangas2/mageckinator/actions/workflows/ci.yml/badge.svg)](https://github.com/ckangas2/mageckinator/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![Container](https://img.shields.io/badge/container-ghcr.io%2Fckangas2%2Fmageckinator-blue)](https://github.com/ckangas2/mageckinator/pkgs/container/mageckinator)

A Nextflow pipeline for MAGeCK reanalysis of a pooled CRISPR screen: library
check → `mageck count` → count QC → `mageck test` (alpha-RRA) → control
calibration → reconciliation against the published hits. Counting and the RRA
test are MAGeCK with standard settings. What this repo adds is the checks
between the steps.

The target dataset is **GSE158298** (Daniloski et al., *Cell* 2021): a
genome-scale GeCKOv2 survival screen in A549-ACE2, SARS-CoV-2 at low (MOI
0.01) and high (0.3) MOI, survivors collected day 6, analysed by MAGeCK RRA
against the plasmid library. It has a published hit list to land near — ACE2,
RAB7A, CTSL, the vacuolar-ATPase and retromer complexes — so the reanalysis
is judged against a fixed target, not explained after the fact.

Each gate below guards a failure that does **not** raise an error: MAGeCK
writes a complete, plausible `gene_summary.txt` whether or not the design, the
counts, or the null are sound. Here, each one stops the run.

## Try it

You need [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 24.04
(which needs Java 17+) and Docker. MAGeCK comes in the container.

```
nextflow run ckangas2/mageckinator -profile test,docker
```

This runs a small synthetic screen (50 targeted genes × 4 guides + 20
nontargeting control pseudo-genes; two Plasmid libraries, three HighMOI
survivor libraries) through the whole pipeline. 8 host-factor genes are
planted enriched ~5× in the survivors, 8 essential genes deplete ~8×, the
rest are flat. `tests/check_test_results.py` then asserts the planted truth
comes back: the host factors as the top positive-selection hits, the
essentials as depleted, the controls null, and the published list recovered.
CI also runs five negative cases, one per gate, and requires each to fail.

No Docker? `-profile test,singularity`, `test,apptainer` or `test,conda`.

## Run it on GSE158298

```
git clone https://github.com/ckangas2/mageckinator && cd mageckinator
bash scripts/fetch_library.sh        # GeCKOv2 human -> reference/GeCKOv2_library.csv
bash scripts/fetch_counts.sh         # fast path: deposited sgRNA count matrix
nextflow run . -profile gse158298,docker \
    --count_table reference/GSE158298_counts.tsv
```

`fetch_counts.sh` pulls the screen's count matrix straight from the GEO
supplementary files, so the reanalysis starts from the same counts the authors
deposited and `mageck count` is skipped. To re-count from reads instead:

```
bash scripts/fetch_runs.sh --list    # how each SRA run's title parses
bash scripts/fetch_runs.sh           # screen FASTQ + samplesheet.csv
nextflow run . -profile gse158298,docker
```

Run accessions are not hand-copied. `fetch_runs.sh` pulls ENA's run report for
the project and `select_runs.py` keeps only the pooled-screen libraries,
assigning each a condition (Plasmid / LowMOI / HighMOI) from its sample title.
Inspect `--list` first; if a title does not parse, fix the pattern, not the
samplesheet.

## Run it on your data

**Samplesheet** (CSV). One row per sequencing library.

```
sample,fastq_1,condition,replicate
Plasmid_rep1,fastq/SRR_plasmid_1.fastq.gz,Plasmid,1
HighMOI_rep1,fastq/SRR_high_1.fastq.gz,HighMOI,1
```

**Comparisons** (CSV). One row per contrast. `selection` is `positive` (hits
enrich in the treatment, e.g. survivors of a killing screen), `negative` (hits
deplete, a dropout screen), or `both`.

```
id,control,treatment,selection
highmoi_vs_plasmid,Plasmid,HighMOI,positive
```

**Library** (CSV): `sgRNA,sequence,gene`, the format `mageck count` reads.
Nontargeting controls should carry a recognisable gene label (`--control_genes`
sets the prefixes; default `Control,NonTargeting,NTC`).

**Run**

```
nextflow run ckangas2/mageckinator -profile docker \
    --input samplesheet.csv --comparisons comparisons.csv --library library.csv \
    --published_hits my_published.tsv   # optional, for reconciliation
```

Already have a count matrix? Skip counting with `--count_table counts.tsv`
(`sgRNA,gene,<samples>`; the sample columns must match `--input`).
`nextflow run ckangas2/mageckinator --help` lists every parameter.

## What gets caught

| failure | what you see | what it actually is | caught by |
|---|---|---|---|
| a comparison names a condition not in the samplesheet | a ranked gene list | MAGeCK was handed an empty/partial group and ranked it anyway | `CHECK_DESIGN` |
| treatment and control are the same condition | "a clean screen, no hits" | there is no contrast; every gene sits at p ≈ 1 | `CHECK_DESIGN` |
| one replicate per arm | confident p-values | RRA variance comes from the sgRNA trend, not the data | `--min_replicates` (`CHECK_DESIGN`) |
| two sgRNAs share a spacer across genes | a full gene ranking | a read maps to one gene arbitrarily; both gene scores are wrong | `VALIDATE_LIBRARY` |
| a gene with one or two guides | a rank like any other | alpha-RRA scores it on one data point | `--min_guides` (`VALIDATE_LIBRARY`) |
| wrong library / un-trimmed adapter | a full count matrix | most reads match nothing and are dropped; DE built on a few percent | mapping-rate floor (`COUNT_QC`) |
| a bottleneck or PCR jackpot | a full count matrix | representation is no longer a fair census | Gini / zero-count / coverage floors (`COUNT_QC`) |
| a plasmid/baseline library swapped or mislabelled | a plausible hit list | the comparison is not the one you think | control-arm replicate correlation (`COUNT_QC`) |
| the null is not calibrated | a long, exciting hit list | normalisation is off; every FDR is optimistic | nontargeting controls must stay out of the hits (`CONTROL_CHECKS`) |
| the selection never took | a clean enrichment list | the screen did not impose the pressure you think | known essentials must deplete (`CONTROL_CHECKS`) |
| treatment/control swapped upstream | a valid ranking | the biology is inverted; `pos` and `neg` columns are switched | opposite direction recovers published hits better (`RECONCILE`) |
| our ranking does not match the paper | a hit list | a different screen, or a pipeline error | recall of the published hits below a floor (`RECONCILE`) |

## Design decisions worth stating

**Control-normalised test.** `mageck test --norm-method control` with the
nontargeting sgRNAs as the normaliser, so a global shift in the targeting
guides (which a killing screen produces) does not drag the null with it. The
control list is fixed once, by `VALIDATE_LIBRARY`, not re-derived per step.

**The null is a check, not decoration.** `CONTROL_CHECKS` requires the
nontargeting guides to behave like a null — roughly uniform two-sided
p-values, essentially none passing FDR. Under `--control-sgrna` MAGeCK uses
those guides as the null and skips them as genes, so the check reads them from
`sgrna_summary.txt`, which scores every control guide, rather than from
`gene_summary.txt`, where they are absent. A screen whose controls are "hits" has a miscalibrated FDR, and no gene
in it can be trusted, however good its rank.

**A positive control, when a list is given.** In a proliferation or survival
screen, core-essential genes drop out. If a curated essential set does not
deplete, the assay did not do what the design assumes, and the enrichment list
is being read off a screen that never worked.

**Selection direction is explicit.** `check_design` records whether hits are
expected to enrich or deplete, and every downstream step reads that rather
than guessing which column of `gene_summary.txt` is the answer. A swapped
contrast inverts the columns silently; `RECONCILE` catches it by testing the
opposite direction against the published hits.

**Gini and zero-count gates apply to the baseline.** Positive-selection
survivors legitimately skew — that is the signal — so those gates default to
the control arms (`--qc_skew_arms`), where representation should still be even.

**Reconciliation is fixed in advance.** The genes to land near
(`assets/daniloski_hits.tsv`) and the recall floor are set before the run, so
agreement is measured, not narrated.

## Prior art, and what this is not

[MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) (Li et al., *Genome
Biol* 2014) does the counting and the RRA test; this pipeline calls it
unchanged. [MAGeCKFlute](https://www.bioconductor.org/packages/MAGeCKFlute/)
does downstream analysis and visualisation, and
[nf-core/crisprseq](https://nf-co.re/crisprseq) is a full, maintained Nextflow
screen pipeline with far more options. For production work, use those.

This repo is narrower: one screen design, end to end, with the checks above
wired in as hard stops, and a validation dataset with a known answer. The part
worth borrowing is the gate layer, not the MAGeCK invocation.

## Outputs

```
results/
  summary.tsv                   per library: reads, mapping, coverage, zero-frac, Gini
  comparisons.tsv               per comparison: hit counts and top genes, in its direction
  design/                       comparisons.validated.tsv, design_report.txt
  library/                      library.validated.csv, control_sgrnas.txt, control_genes.txt
  count/                        counts.count.txt, counts.countsummary.txt, count_qc.tsv
  test/<id>/                    <id>.gene_summary.txt, <id>.sgrna_summary.txt
  qc/                           <id>.control_checks.tsv, <id>.reconciliation.tsv
  pipeline_info/                Nextflow report, timeline, trace, DAG
```

## Layout

```
main.nf                    design -> library -> count -> QC gate -> test -> controls -> reconcile
nextflow.config            parameters, resources, profiles
conf/                      test.config, gse158298.config
modules/local/             one process per file, each with a stub block
bin/
  check_design.py          samplesheet + comparisons -> validated design
  validate_library.py      sgRNA library sanity; control list
  table_to_countsummary.py countsummary from a deposited count matrix
  count_qc.py              mapping / coverage / Gini / baseline-swap gate
  control_checks.py        nontargeting null + essential-gene dropout
  reconcile.py             our ranking vs the published hits
  collect_summary.py       run-level tables
scripts/                   library, count-matrix and FASTQ download; run selection
assets/                    daniloski_hits.tsv, essential_genes_core.txt, comparisons
tests/
  make_test_data.py        regenerates tests/data/ byte-for-byte (CI checks this)
  test_checks.py           unit tests: one planted failure per gate
  check_test_results.py    asserts the test profile against the planted truth
```

## Development

```
bash tests/run_tests.sh                          # unit tests, seconds, no MAGeCK
nextflow run . -profile test,docker -stub-run    # wiring only, no MAGeCK
nextflow run . -profile test,docker              # end to end
python3 tests/check_test_results.py results
```

`-stub-run` fabricates the MAGeCK outputs but still runs the Python gates, so
outside a container it needs `numpy` (`pip install numpy pandas`). The unit
tests need `numpy pandas pytest` and nothing else.

## Status

- [x] Unit tests (20) pass; `nextflow lint` clean on 26.10.4; stub run completes, and design/library failures stop it before counting
- [x] `-profile test,docker` end to end locally (Nextflow 26.10.4): real MAGeCK recovers all 8 planted host factors (pos FDR<0.05), all 8 essentials depleted, controls null, published recall 1.00
- [ ] The same in CI, plus the five negative cases
- [ ] Container build and conda solve in CI
- [ ] GSE158298 run; RRA hits vs `assets/daniloski_hits.tsv`; control calibration; reconciliation → RESULTS.md
