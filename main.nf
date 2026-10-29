#!/usr/bin/env nextflow
/*
 * mageckinator — MAGeCK reanalysis of a pooled CRISPR screen, with the checks
 * that stand between a count table and a trustworthy hit list.
 *
 *   design:   samplesheet + comparisons validated before any counting
 *   library:  sgRNA library sanity; control-sgRNA list fixed here
 *   counts:   mageck count (FASTQ) OR a deposited count table
 *   gate:     count QC (mapping, coverage, Gini, baseline-swap detection)
 *   test:     mageck test (alpha-RRA), control-normalised, per comparison
 *   controls: nontargeting null calibrated; essentials deplete
 *   reconcile: our ranking vs the published hits (optional)
 *
 * Every gate stops the run. Each guards a failure that otherwise produces a
 * complete, plausible gene_summary.txt with the wrong genes in it.
 */
nextflow.enable.dsl = 2

include { CHECK_DESIGN      } from './modules/local/check_design'
include { VALIDATE_LIBRARY  } from './modules/local/validate_library'
include { MAGECK_COUNT      } from './modules/local/mageck_count'
include { COUNTS_FROM_TABLE } from './modules/local/counts_from_table'
include { COUNT_QC          } from './modules/local/count_qc'
include { MAGECK_TEST       } from './modules/local/mageck_test'
include { CONTROL_CHECKS    } from './modules/local/control_checks'
include { RECONCILE         } from './modules/local/reconcile'
include { COLLECT_SUMMARY   } from './modules/local/collect_summary'

def helpMessage() {
    log.info """
    mageckinator v${workflow.manifest.version}

    Usage:
      nextflow run ckangas2/mageckinator -profile test,docker
      nextflow run ckangas2/mageckinator -profile docker \\
          --input samplesheet.csv --comparisons comparisons.csv \\
          --library library.csv

    Input
      --input          CSV: sample,fastq_1,condition,replicate
      --comparisons    CSV: id,control,treatment,selection
                       (selection: positive | negative | both)
      --library        sgRNA library CSV: sgRNA,sequence,gene

    Counts (choose one)
      (default)        count FASTQ in --input with `mageck count`
      --count_table    start from a deposited matrix: sgRNA,gene,<samples>
                       (sample columns must match --input sample names)

    Gates
      --min_replicates / --min_control_replicates
      --min_guides             sgRNAs per targeted gene (default ${params.min_guides})
      --min_mapped_frac / --min_coverage / --max_zero_frac / --max_gini
      --min_replicate_rho      control-arm replicate correlation floor
      --max_control_hit_frac / --min_control_median_p
      --essential_genes        file of essential ids (positive control; optional)
      --min_essential_hit_frac
      --published_hits         gene[,rank] list to reconcile against (optional)
      --reconcile_min_recall / --reconcile_top_n

    MAGeCK
      --control_genes          gene-label prefixes for controls (default '${params.control_genes}')
      --sgrna_len              spacer length for mageck count (default ${params.sgrna_len})
      --fdr                    gene FDR for hit counts (default ${params.fdr})

    Output
      --outdir (default: ${params.outdir})
    """.stripIndent()
}

def isSet(v) {
    return v != null && !(v instanceof Boolean) && !(v.toString().trim() in ['', 'null', 'false', 'true'])
}

def validateParams() {
    def problems = []
    if (!isSet(params.input))       problems << "--input is required"
    if (!isSet(params.comparisons)) problems << "--comparisons is required"
    if (!isSet(params.library))     problems << "--library is required"
    if (problems) {
        problems.each { log.error it }
        error "Fix the parameters above (see --help)."
    }
}

def resolvePath(base, String p) {
    return (p.contains('://') || p.startsWith('/')) ? p : "${base}/${p}"
}

workflow {
    if (params.help) { helpMessage(); return }
    validateParams()
    main_flow()
}

workflow main_flow {
    main:
    def samplesheet = file(params.input,       checkIfExists: true)
    def comparisons = file(params.comparisons, checkIfExists: true)
    def library     = file(params.library,     checkIfExists: true)
    def base        = samplesheet.parent

    def rows = samplesheet.splitCsv(header: true, strip: true)
    if (rows.size() == 0) error "samplesheet ${params.input} has a header but no samples"

    // ---- design and library, both before any counting
    CHECK_DESIGN(samplesheet, comparisons)
    VALIDATE_LIBRARY(library)

    // ---- counts: from a deposited table, or from FASTQ
    def counts
    def countsummary
    if (isSet(params.count_table)) {
        COUNTS_FROM_TABLE(file(params.count_table, checkIfExists: true))
        counts       = COUNTS_FROM_TABLE.out.counts
        countsummary = COUNTS_FROM_TABLE.out.summary
    } else {
        def reads = channel.fromList(rows).map { row ->
            tuple(row.sample, file(resolvePath(base, row.fastq_1), checkIfExists: true))
        }
        def sorted       = reads.toSortedList { a, b -> a[0] <=> b[0] }
        def count_files  = sorted.map { lst -> lst.collect { it[1] } }
        def count_labels = sorted.map { lst -> lst.collect { it[0] }.join(',') }
        MAGECK_COUNT(count_files, count_labels,
                     VALIDATE_LIBRARY.out.library, VALIDATE_LIBRARY.out.control_sgrnas)
        counts       = MAGECK_COUNT.out.counts
        countsummary = MAGECK_COUNT.out.summary
    }

    // ---- count QC gate
    COUNT_QC(countsummary, counts, samplesheet, CHECK_DESIGN.out.comparisons)

    // ---- per comparison, gated on COUNT_QC passing
    def gate = COUNT_QC.out.qc.map { true }
    def comparison_ch = CHECK_DESIGN.out.comparisons
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.id, row.control_samples, row.treatment_samples, row.selection) }
        .combine(gate)
        .map { id, c, t, s, _g -> tuple(id, c, t, s) }

    // counts/controls are value channels (value-channel inputs upstream), so
    // they are reusable across every comparison without .first().
    MAGECK_TEST(comparison_ch, counts, VALIDATE_LIBRARY.out.control_sgrnas)

    // ---- control calibration (every comparison)
    def essential = isSet(params.essential_genes) ?
        file(params.essential_genes, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    CONTROL_CHECKS(MAGECK_TEST.out.gene_summary,
                   VALIDATE_LIBRARY.out.control_genes, essential)

    // ---- reconcile against the published hits (optional)
    if (isSet(params.published_hits)) {
        RECONCILE(MAGECK_TEST.out.gene_summary,
                  file(params.published_hits, checkIfExists: true))
    }

    // ---- run-level tables
    def gs       = MAGECK_TEST.out.gene_summary
    def gs_files = gs.map { id, sel, f -> f }.collect()
    def gs_args  = gs.map { id, sel, f -> "${id}=${f.name}" }.collect()
    COLLECT_SUMMARY(COUNT_QC.out.qc, CHECK_DESIGN.out.comparisons, gs_files, gs_args)

    COLLECT_SUMMARY.out.comparisons.subscribe { _f ->
        log.info "mageckinator: per-library QC in ${params.outdir}/summary.tsv, " +
                 "per-comparison hits in ${params.outdir}/comparisons.tsv"
    }
}
