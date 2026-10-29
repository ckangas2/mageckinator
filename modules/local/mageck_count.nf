// FASTQ -> sgRNA count matrix with `mageck count`. The library header line is
// stripped first: mageck count treats every line as a record, so a header
// would be loaded as a junk sgRNA named "sgRNA" with sequence "sequence".
process MAGECK_COUNT {
    tag "count"
    label 'process_medium'
    publishDir "${params.outdir}/count", mode: params.publish_dir_mode,
        pattern: '*.{txt,log,count_normalized.txt}'

    input:
    path fastqs, stageAs: 'fq/*'
    val  labels
    path library
    path control_sgrnas

    output:
    path 'counts.count.txt',        emit: counts
    path 'counts.countsummary.txt', emit: summary
    path 'counts.log',              emit: log, optional: true

    script:
    def files = fastqs instanceof List ? fastqs.join(' ') : "${fastqs}"
    def ctrl = params.sgrna_len ? "--sgrna-len ${params.sgrna_len}" : ''
    """
    tail -n +2 ${library} > library.noheader.csv
    mageck count \\
        --list-seq library.noheader.csv \\
        --fastq ${files} \\
        --sample-label ${labels} \\
        --output-prefix counts ${ctrl} \\
        --control-sgrna ${control_sgrnas} 2> counts.log
    """

    stub:
    """
    # Minimal count table + summary so -stub-run exercises the wiring without mageck.
    python3 - "${labels}" ${library} <<'PY'
import sys
labels = sys.argv[1].split(',')
lib = [l.strip().split(',') for l in open(sys.argv[2]) if l.strip()][1:]
with open('counts.count.txt','w') as f:
    f.write('sgRNA\\tGene\\t'+'\\t'.join(labels)+'\\n')
    for sg,seq,gene in lib:
        f.write(sg+'\\t'+gene+'\\t'+'\\t'.join('100' for _ in labels)+'\\n')
with open('counts.countsummary.txt','w') as f:
    f.write('File\\tLabel\\tReads\\tMapped\\tPercentage\\tTotalsgRNAs\\tZerocounts\\tGiniIndex\\n')
    for s in labels:
        f.write(f'{s}.fq\\t{s}\\t{len(lib)*100}\\t{len(lib)*100}\\t0.95\\t{len(lib)}\\t0\\t0.10\\n')
PY
    touch counts.log
    """
}
