# Results: GSE158298

Not run yet. This file is the template the GSE158298 run fills in, with the
comparison points chosen in advance so the result is judged against them, not
explained after.

Daniloski et al., *Cell* 2021, "Identification of Required Host Factors for
SARS-CoV-2 Infection in Human Cells." A549-ACE2, GeCKOv2 (122,411 sgRNAs,
19,050 genes, ~1,000 nontargeting controls), SARS-CoV-2 USA-WA1/2020 at MOI
0.01 and 0.3, survivors collected day 6, MAGeCK RRA vs the plasmid library.
GEO: GSE158298.

## What to compare against

**Hits (positive selection).** The authors report ~1,000 genes with
significant RRA enrichment (p < 0.05) and a smaller high-confidence set
assigned to the viral life cycle. The genes to land near are in
`assets/daniloski_hits.tsv`:

| stage | genes |
|---|---|
| attachment / entry | ACE2, CTSL |
| actin / ARP2-3 | ACTR2, ACTR3, ARPC3, ARPC4 |
| vacuolar-ATPase | ATP6AP1, ATP6AP2, ATP6V0B/0C/0D1, ATP6V1A/1B2/1C1/1E1/1G1/1H, TMEM199, TOR1AIP1 |
| endosome / retromer | RAB7A, VPS35, VPS29, VPS26A, CCDC22 |
| other | NPC1, PIK3C3 |

ACE2 should be at or near the top of the positive-selection list in both MOI
arms. The paper notes substantial overlap of enriched genes between the low
and high MOI conditions — report that overlap.

**Reconciliation.** `reconcile.py` reports recall of this list in our top-N
and below our FDR, for each comparison. Recall will be below 1.0 — a different
MAGeCK version and count-vs-FASTQ starting point move ranks — but the vacuolar-
ATPase and retromer members should come back as a block, and ACE2/RAB7A/CTSL
near the top. If the opposite (negative) direction recovers them better, the
contrast is inverted; stop and check the design before interpreting.

**Controls.** Nontargeting controls must stay out of the hit list and sit near
a uniform p-value distribution (`control_checks.tsv`). Core-essential genes
(`assets/essential_genes_core.txt`, or Hart CEGv2 via `--essential_genes`)
should deplete — the screen's positive control.

**Count QC.** The paper reports high plasmid↔post-selection correlation
(r ≈ 0.84) and much greater guide dropout after infection (expected: survivors
are a bottleneck). So expect the baseline arms to pass the Gini/zero gates and
the survivor arms to skew — which is why those gates default to the baseline.

## Sections to fill

1. Run QC: `summary.tsv` — reads, mapping rate, coverage, zero-count fraction,
   Gini per library; plasmid-replicate correlation.
2. Design report: samples per comparison, conditions, what was excluded.
3. Hits: counts at FDR < 0.05 and at p < 0.05 for each MOI arm, top genes,
   low-vs-high MOI overlap. (The paper thresholded on p; report both.)
4. Control calibration: nontargeting null (control hit fraction, median
   p-value), essential-gene depletion.
5. Reconciliation table vs `assets/daniloski_hits.tsv`, per comparison.
6. Fast path (deposited counts) vs full path (re-counted FASTQ): where the two
   rankings agree and where they diverge.
7. Anything that disagreed with expectations, and what it turned out to be.
