// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process MERGE_COUNTS_RSEM {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path gene_counts
    path isoform_counts
    val  options

    output:
    path "rsem.merged.gene_tpm.tsv"         , emit: tpm_gene
    path "rsem.merged.gene_counts.tsv"      , emit: counts_gene
    path "rsem.merged.transcript_tpm.tsv"   , emit: tpm_transcript
    path "rsem.merged.transcript_counts.tsv", emit: counts_transcript

    script:
    """
    mkdir tmp_genes
    echo "gene_id\tgene_symbol" > gene_ids.txt
    cut -f 1 ${gene_counts[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> gene_ids.txt
    for fileid in $gene_counts; do
        basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.tpm.txt
        grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_genes/\${fileid}.tpm.txt
        basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.counts.txt
        grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_genes/\${fileid}.counts.txt
    done

    mkdir tmp_isoforms
    echo "transcript_id\tgene_symbol" > transcript_ids.txt
    cut -f 1 ${isoform_counts[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> transcript_ids.txt
    for fileid in $isoform_counts; do
        basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.tpm.txt
        grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_isoforms/\${fileid}.tpm.txt
        basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.counts.txt
        grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_isoforms/\${fileid}.counts.txt
    done

    paste gene_ids.txt tmp_genes/*.tpm.txt > rsem.merged.gene_tpm.tsv
    paste gene_ids.txt tmp_genes/*.counts.txt > rsem.merged.gene_counts.tsv
    paste transcript_ids.txt tmp_isoforms/*.tpm.txt > rsem.merged.transcript_tpm.tsv
    paste transcript_ids.txt tmp_isoforms/*.counts.txt > rsem.merged.transcript_counts.tsv
    """
}
