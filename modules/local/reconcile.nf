// Our ranking vs the genes the original study reported. Recall below a floor,
// or the opposite selection recovering the published hits better (an inverted
// comparison), stops the run.
process RECONCILE {
    tag "${id}"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    tuple val(id), val(selection), path(gene_summary)
    path published

    output:
    tuple val(id), path("${id}.reconciliation.tsv"), emit: report

    script:
    def sel = selection == 'negative' ? 'negative' : 'positive'
    """
    reconcile.py \\
        --gene-summary ${gene_summary} \\
        --published ${published} \\
        --selection ${sel} \\
        --fdr ${params.fdr} \\
        --top-n ${params.reconcile_top_n} \\
        --min-recall ${params.reconcile_min_recall} \\
        --out ${id}.reconciliation.tsv
    """

    stub:
    """
    reconcile.py --gene-summary ${gene_summary} --published ${published} \\
        --selection ${selection == 'negative' ? 'negative' : 'positive'} \\
        --min-recall 0.0 --out ${id}.reconciliation.tsv
    """
}
