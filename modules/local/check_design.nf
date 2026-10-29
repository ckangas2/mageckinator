// Validate samplesheet + comparisons into a frozen design, before any counting.
// A bad comparison is caught here in milliseconds, not after a genome-scale
// count has run.
process CHECK_DESIGN {
    tag "design"
    label 'process_single'
    publishDir "${params.outdir}/design", mode: params.publish_dir_mode

    input:
    path samplesheet, stageAs: 'in/*'
    path comparisons, stageAs: 'in/*'

    output:
    path 'comparisons.validated.tsv', emit: comparisons
    path 'design_report.txt',         emit: report

    script:
    def ctrl = params.min_control_replicates ? "--min-control-replicates ${params.min_control_replicates}" : ''
    """
    check_design.py \\
        --samplesheet ${samplesheet} \\
        --comparisons ${comparisons} \\
        --min-replicates ${params.min_replicates} ${ctrl} \\
        --out-dir .
    """

    stub:
    """
    check_design.py --samplesheet ${samplesheet} --comparisons ${comparisons} \\
        --min-replicates ${params.min_replicates} --out-dir .
    """
}
