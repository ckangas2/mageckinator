// Did the screen behave? Nontargeting controls must look like the null, and
// (if an essential-gene list is given) known essentials must deplete. Either
// failure stops the comparison: the hit list cannot be trusted.
process CONTROL_CHECKS {
    tag "${id}"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    tuple val(id), val(selection), path(gene_summary)
    path control_genes
    path essential_genes

    output:
    tuple val(id), path("${id}.control_checks.tsv"), emit: report

    script:
    def ess = essential_genes.name != 'NO_FILE' ? "--essential-genes ${essential_genes}" : ''
    """
    control_checks.py \\
        --gene-summary ${gene_summary} \\
        --control-genes '${params.control_genes}' \\
        --control-genes-list ${control_genes} \\
        --fdr ${params.fdr} \\
        --max-control-hit-frac ${params.max_control_hit_frac} \\
        --min-control-median-p ${params.min_control_median_p} ${ess} \\
        --min-essential-hit-frac ${params.min_essential_hit_frac} \\
        --out ${id}.control_checks.tsv
    """

    stub:
    """
    control_checks.py --gene-summary ${gene_summary} \\
        --control-genes-list ${control_genes} \\
        --max-control-hit-frac 1.0 --min-control-median-p 0.0 \\
        --out ${id}.control_checks.tsv
    """
}
