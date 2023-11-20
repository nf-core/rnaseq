process SUMMARIZEDEXPERIMENT {
    tag "$tx2gene"
    label "process_medium"

    conda "bioconda::bioconductor-summarizedexperiment=1.24.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.24.0--r41hdfd78af_0' :
        'biocontainers/bioconductor-summarizedexperiment:1.24.0--r41hdfd78af_0' }"

    input:
    path counts
    path tpm
    path tx2gene

    output:
    path "*.rds"       , emit: rds
    tuple val("${task.process}"), val('r-base'), cmd("echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//'"), emit: versions1
    tuple val("${task.process}"), val('bioconductor-summarizedexperiment'), cmd("Rscript -e 'library(SummarizedExperiment); cat(as.character(packageVersion(\'SummarizedExperiment\')))'"), emit: versions2

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    summarizedexperiment.r \\
        NULL \\
        $counts \\
        $tpm \\
        $tx2gene
    """
}
