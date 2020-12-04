// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SALMON_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "biocontainers/biocontainers:v1.2.0_cv1"
    
    input:
    path ('genes_counts/*')
    path ('genes_tpm/*')
    path ('genes_counts_length_scaled/*')
    path ('genes_tpm_length_scaled/*')
    path ('genes_counts_scaled/*')
    path ('genes_tpm_scaled/*')
    path ('isoforms_counts/*')
    path ('isoforms_tpm/*')
    
    output:
    path "salmon.merged.gene_counts.tsv"              , emit: counts_gene
    path "salmon.merged.gene_tpm.tsv"                 , emit: tpm_gene
    path "salmon.merged.gene_counts_length_scaled.tsv", emit: counts_gene_length_scaled
    path "salmon.merged.gene_tpm_length_scaled.tsv"   , emit: tpm_gene_length_scaled
    path "salmon.merged.gene_counts_scaled.tsv"       , emit: counts_gene_scaled
    path "salmon.merged.gene_tpm_scaled.tsv"          , emit: tpm_gene_scaled
    path "salmon.merged.transcript_counts.tsv"        , emit: counts_transcript
    path "salmon.merged.transcript_tpm.tsv"           , emit: tpm_transcript

    script:
    """
    mkdir -p tmp/genes_counts
    echo "${params.fc_group_features}" > gene_ids.txt
    cut -f 1 `ls ./genes_counts/* | head -n 1` | tail -n +2 >> gene_ids.txt
    for fileid in `ls ./genes_counts/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_counts/\${filename}
    done

    mkdir -p tmp/genes_tpm
    for fileid in `ls ./genes_tpm/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_tpm/\${filename}
    done

    mkdir -p tmp/genes_counts_length_scaled
    for fileid in `ls ./genes_counts_length_scaled/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_counts_length_scaled/\${filename}
    done

    mkdir -p tmp/genes_tpm_length_scaled
    for fileid in `ls ./genes_tpm_length_scaled/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_tpm_length_scaled/\${filename}
    done

    mkdir -p tmp/genes_counts_scaled
    for fileid in `ls ./genes_counts_scaled/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_counts_scaled/\${filename}
    done

    mkdir -p tmp/genes_tpm_scaled
    for fileid in `ls ./genes_tpm_scaled/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/genes_tpm_scaled/\${filename}
    done

    mkdir -p tmp/isoforms_counts
    echo "transcript_id" > transcript_ids.txt
    cut -f 1 `ls ./isoforms_counts/* | head -n 1` | tail -n +2 >> transcript_ids.txt
    for fileid in `ls ./isoforms_counts/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/isoforms_counts/\${filename}
    done

    mkdir -p tmp/isoforms_tpm
    for fileid in `ls ./isoforms_tpm/*`; do
        filename=`basename \$fileid`
        cut -f 2 \${fileid} > tmp/isoforms_tpm/\${filename}
    done

    paste gene_ids.txt tmp/genes_counts/* > salmon.merged.gene_counts.tsv
    paste gene_ids.txt tmp/genes_tpm/* > salmon.merged.gene_tpm.tsv
    paste gene_ids.txt tmp/genes_counts_length_scaled/* > salmon.merged.gene_counts_length_scaled.tsv
    paste gene_ids.txt tmp/genes_tpm_length_scaled/* > salmon.merged.gene_tpm_length_scaled.tsv
    paste gene_ids.txt tmp/genes_counts_scaled/* > salmon.merged.gene_counts_scaled.tsv
    paste gene_ids.txt tmp/genes_tpm_scaled/* > salmon.merged.gene_tpm_scaled.tsv
    paste transcript_ids.txt tmp/isoforms_counts/* > salmon.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms_tpm/* > salmon.merged.transcript_tpm.tsv
    """
}
