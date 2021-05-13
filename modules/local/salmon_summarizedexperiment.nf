// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SALMON_SUMMARIZEDEXPERIMENT {
    tag "$tx2gene"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-summarizedexperiment=1.20.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.20.0--r40_0"
    } else {
        container "quay.io/biocontainers/bioconductor-summarizedexperiment:1.20.0--r40_0"
    }

    input:
    path counts
    path tpm
    path tx2gene

    output:
    path "*.rds"         , emit: rds
    path  "*.version.txt", emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_summarizedexperiment.r NULL $counts $tpm
    Rscript -e "library(SummarizedExperiment); write(x=as.character(packageVersion('SummarizedExperiment')), file='bioconductor-summarizedexperiment.version.txt')"
    """
}
