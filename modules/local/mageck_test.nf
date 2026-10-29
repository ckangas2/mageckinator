// One comparison through `mageck test` (alpha-RRA), normalised on the
// nontargeting controls so the null and the FDR are control-calibrated.
process MAGECK_TEST {
    tag "${id}"
    label 'process_low'
    publishDir "${params.outdir}/test/${id}", mode: params.publish_dir_mode

    input:
    tuple val(id), val(control_label), val(treatment_label), val(selection)
    path counts
    path control_sgrnas

    output:
    tuple val(id), val(selection), path("${id}.gene_summary.txt"), emit: gene_summary
    path "${id}.sgrna_summary.txt",                                emit: sgrna_summary
    path "${id}.log",                                              emit: log, optional: true

    script:
    """
    mageck test \\
        --count-table ${counts} \\
        --treatment-id ${treatment_label} \\
        --control-id ${control_label} \\
        --control-sgrna ${control_sgrnas} \\
        --norm-method control \\
        --output-prefix ${id} 2> ${id}.log
    """

    stub:
    """
    python3 - ${id} ${counts} <<'PY'
import sys
cid, counts = sys.argv[1], sys.argv[2]
genes, seen = [], set()
for i, line in enumerate(open(counts)):
    if i == 0:
        continue
    g = line.split('\\t')[1]
    if g not in seen:
        seen.add(g); genes.append(g)
cols = ['id','num','neg|score','neg|p-value','neg|fdr','neg|rank','neg|goodsgrna','neg|lfc',
        'pos|score','pos|p-value','pos|fdr','pos|rank','pos|goodsgrna','pos|lfc']
with open(f'{cid}.gene_summary.txt','w') as f:
    f.write('\\t'.join(cols)+'\\n')
    for r,g in enumerate(genes,1):
        f.write('\\t'.join([g,'4','0.5','0.5','0.9',str(r),'0','0','0.5','0.5','0.9',str(r),'0','0'])+'\\n')
open(f'{cid}.sgrna_summary.txt','w').write('sgrna\\tGene\\n')
PY
    touch ${id}.log
    """
}
