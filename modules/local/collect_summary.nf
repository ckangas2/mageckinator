// Assemble the two run-level tables a reader opens first: summary.tsv (per
// library) and comparisons.tsv (hits per comparison, in its own direction).
process COLLECT_SUMMARY {
    tag "summary"
    label 'process_single'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path count_qc
    path validated
    path gene_summaries, stageAs: 'gs/*'
    val  gene_summary_args

    output:
    path 'summary.tsv',      emit: summary
    path 'comparisons.tsv',  emit: comparisons

    script:
    def args = gene_summary_args.collect { it.replaceFirst('=', '=gs/') }.join(' ')
    """
    collect_summary.py \\
        --count-qc ${count_qc} \\
        --comparisons ${validated} \\
        --gene-summary ${args} \\
        --fdr ${params.fdr} \\
        --top ${params.summary_top} \\
        --out-dir .
    """

    stub:
    """
    collect_summary.py --count-qc ${count_qc} --comparisons ${validated} \\
        --gene-summary ${gene_summary_args.collect { it.replaceFirst('=', '=gs/') }.join(' ')} \\
        --out-dir .
    """
}
