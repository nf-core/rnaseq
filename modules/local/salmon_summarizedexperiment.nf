process SALMON_SUMMARIZEDEXPERIMENT {
    tag "$tx2gene"
    label "process_medium"

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
    path "*.rds"       , emit: rds
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_summarizedexperiment.r \\
        NULL \\
        $counts \\
        $tpm

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
