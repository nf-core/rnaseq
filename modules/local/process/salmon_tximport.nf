// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SALMON_TXIMPORT {
    //tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bioconductor-tximeta=1.8.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0"
    } else {
        container "quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0"
    }

    input:
    path ("salmon/*")
    path  tx2gene
    
    output:
    path("*gene_tpm.tsv")                 , emit: tpm_gene
    path("*gene_counts.tsv")              , emit: counts_gene
    path("*gene_tpm_length_scaled.tsv")   , emit: tpm_gene_length_scaled
    path("*gene_counts_length_scaled.tsv"), emit: counts_gene_length_scaled
    path("*gene_tpm_scaled.tsv")          , emit: tpm_gene_scaled
    path("*gene_counts_scaled.tsv")       , emit: counts_gene_scaled
    path("*transcript_tpm.tsv")           , emit: tpm_transcript
    path("*transcript_counts.tsv")        , emit: counts_transcript
    path  "*.version.txt"                                  , emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tximport.r NULL salmon salmon.merged
    Rscript -e "library(tximeta); write(x=as.character(packageVersion('tximeta')), file='bioconductor-tximeta.version.txt')"
    """
}
