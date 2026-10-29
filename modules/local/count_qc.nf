// Gate on the count table: mapping rate, coverage, zero fraction, Gini, and
// control-arm replicate correlation (swap detection). Stops the run before any
// gene-level test if a library is not what it claims to be.
process COUNT_QC {
    tag "count-qc"
    label 'process_single'
    publishDir "${params.outdir}/count", mode: params.publish_dir_mode

    input:
    path countsummary
    path counts
    path samplesheet, stageAs: 'in/*'
    path validated

    output:
    path 'count_qc.tsv', emit: qc

    script:
    """
    count_qc.py \\
        --countsummary ${countsummary} \\
        --counts ${counts} \\
        --samplesheet ${samplesheet} \\
        --comparisons ${validated} \\
        --min-mapped-frac ${params.min_mapped_frac} \\
        --max-zero-frac ${params.max_zero_frac} \\
        --max-gini ${params.max_gini} \\
        --min-coverage ${params.min_coverage} \\
        --min-replicate-rho ${params.min_replicate_rho} \\
        --skew-arms ${params.qc_skew_arms} \\
        --out count_qc.tsv
    """

    stub:
    """
    # wiring only: disable the thresholds so -stub-run always completes
    count_qc.py --countsummary ${countsummary} --counts ${counts} \\
        --samplesheet ${samplesheet} --comparisons ${validated} \\
        --min-coverage 1 --min-mapped-frac 0 --min-replicate-rho 0 \\
        --max-gini 1 --max-zero-frac 1 --out count_qc.tsv
    """
}
