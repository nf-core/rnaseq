// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process SALMON_MERGE_COUNTS {
    tag "$tx2gene"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path gene_tpm
    path gene_counts
    path transcript_tpm
    path transcript_counts
    path tx2gene
    val  options

    output:
    path "salmon.merged.gene_tpm.csv"         , emit: tpm_gene
    path "salmon.merged.gene_counts.csv"      , emit: counts_gene
    path "salmon.merged.transcript_tpm.csv"   , emit: tpm_transcript
    path "salmon.merged.transcript_counts.csv", emit: counts_transcript
    path "*.rds"                              , emit: rds

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    // First field is the gene/transcript ID
    gene_ids       = "<(cut -f1 -d, ${gene_tpm[0]} | tail -n +2 | cat <(echo '${params.fc_group_features}') - )"
    transcript_ids = "<(cut -f1 -d, ${transcript_tpm[0]} | tail -n +2 | cat <(echo 'transcript_id') - )"

    // Second field is counts/TPM
    gene_tpm_cols          = gene_tpm.collect { f -> "<(cut -d, -f2 ${f})" }.join(" ")
    gene_counts_cols       = gene_counts.collect { f -> "<(cut -d, -f2 ${f})" }.join(" ")
    transcript_tpm_cols    = transcript_tpm.collect { f -> "<(cut -d, -f2 ${f})" }.join(" ")
    transcript_counts_cols = transcript_counts.collect { f -> "<(cut -d, -f2 ${f})" }.join(" ")
    """
    paste -d, $gene_ids $gene_tpm_cols > salmon.merged.gene_tpm.csv
    paste -d, $gene_ids $gene_counts_cols > salmon.merged.gene_counts.csv
    paste -d, $transcript_ids $transcript_tpm_cols > salmon.merged.transcript_tpm.csv
    paste -d, $transcript_ids $transcript_counts_cols > salmon.merged.transcript_counts.csv

    salmon_se.r NULL salmon.merged.gene_counts.csv salmon.merged.gene_tpm.csv
    salmon_se.r NULL salmon.merged.transcript_counts.csv salmon.merged.transcript_tpm.csv
    """
}
