// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process RSEM_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "biocontainers/biocontainers:v1.2.0_cv1"
    
    conda (params.conda ? "conda-forge::sed=4.7" : null)
    
    input:
    path ('genes/*')
    path ('isoforms/*')
    val  options

    output:
    path "rsem.merged.gene_counts.tsv"      , emit: counts_gene
    path "rsem.merged.gene_tpm.tsv"         , emit: tpm_gene
    path "rsem.merged.transcript_counts.tsv", emit: counts_transcript
    path "rsem.merged.transcript_tpm.tsv"   , emit: tpm_transcript
    
    script:
    """
    mkdir tmp_genes
    echo "gene_id\tgene_symbol" > gene_ids.txt
    cut -f 1 `ls ./genes/* | head -n 1` | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> gene_ids.txt
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        echo \$samplename > tmp_genes/\${samplename}.counts.txt
        grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_genes/\${samplename}.counts.txt
        echo \$samplename > tmp_genes/\${samplename}.tpm.txt
        grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_genes/\${samplename}.tpm.txt
    done

    mkdir tmp_isoforms
    echo "transcript_id\tgene_symbol" > transcript_ids.txt
    cut -f 1 `ls ./isoforms/* | head -n 1` | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> transcript_ids.txt
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        echo \$samplename > tmp_isoforms/\${samplename}.counts.txt
        grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_isoforms/\${samplename}.counts.txt
        echo \$samplename > tmp_isoforms/\${samplename}.tpm.txt
        grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_isoforms/\${samplename}.tpm.txt
    done

    paste gene_ids.txt tmp_genes/*.counts.txt > rsem.merged.gene_counts.tsv
    paste gene_ids.txt tmp_genes/*.tpm.txt > rsem.merged.gene_tpm.tsv
    paste transcript_ids.txt tmp_isoforms/*.counts.txt > rsem.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp_isoforms/*.tpm.txt > rsem.merged.transcript_tpm.tsv
    """
}
