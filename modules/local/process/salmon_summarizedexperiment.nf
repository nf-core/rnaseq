// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process SALMON_SUMMARIZEDEXPERIMENT {
    tag "$tx2gene"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/bioconductor-summarizedexperiment:1.18.1--r40_0"
    //container https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.18.1--r40_0

    conda (params.conda ? "bioconda::bioconductor-summarizedexperiment=1.18.1" : null)

    input:
    path counts
    path tpm
    path tx2gene
    val  options

    output:
    path "*.rds"         , emit: rds
    path  "*.version.txt", emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_summarizedexperiment.r NULL $counts $tpm
    Rscript -e "library(SummarizedExperiment); write(x=as.character(packageVersion('SummarizedExperiment')), file='bioconductor-summarizedexperiment.version.txt')"
    """
}
