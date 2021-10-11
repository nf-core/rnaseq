process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::samtools=1.13' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
    } else {
        container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path  "versions.yml"               , emit: versions

    script:
    """
    samtools flagstat $bam > ${bam}.flagstat
    cat <<-END_VERSIONS > versions.yml
    SAMTOOLS_FLAGSTAT:
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
