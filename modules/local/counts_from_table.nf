// Start from a deposited count matrix instead of FASTQ. The table becomes the
// count file; a countsummary is derived from it (mapping rate unknown, so the
// mapping gate is a no-op on this path -- see table_to_countsummary.py).
process COUNTS_FROM_TABLE {
    tag "count-table"
    label 'process_single'
    publishDir "${params.outdir}/count", mode: params.publish_dir_mode

    input:
    path table, stageAs: 'in/*'

    output:
    path 'counts.count.txt',        emit: counts
    path 'counts.countsummary.txt', emit: summary

    script:
    """
    cp ${table} counts.count.txt
    table_to_countsummary.py --counts counts.count.txt --out counts.countsummary.txt
    """

    stub:
    """
    cp ${table} counts.count.txt
    table_to_countsummary.py --counts counts.count.txt --out counts.countsummary.txt
    """
}
