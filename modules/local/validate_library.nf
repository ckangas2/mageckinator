// sgRNA library sanity, before counting. Writes a cleaned, MAGeCK-ready
// library and the authoritative control-sgRNA list used for normalisation and
// the null.
process VALIDATE_LIBRARY {
    tag "library"
    label 'process_single'
    publishDir "${params.outdir}/library", mode: params.publish_dir_mode

    input:
    path library, stageAs: 'in/*'

    output:
    path 'library.validated.csv', emit: library
    path 'control_sgrnas.txt',    emit: control_sgrnas
    path 'control_genes.txt',     emit: control_genes
    path 'library_report.txt',    emit: report

    script:
    """
    validate_library.py \\
        --library ${library} \\
        --control-genes '${params.control_genes}' \\
        --min-guides ${params.min_guides} \\
        --out-dir .
    """

    stub:
    """
    validate_library.py --library ${library} \\
        --control-genes '${params.control_genes}' \\
        --min-guides ${params.min_guides} --out-dir .
    """
}
