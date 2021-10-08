// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0"
    } else {
        container "quay.io/biocontainers/gffread:0.12.1--h8b12597_0"
    }

    input:
    path gff

    output:
    path "*.gtf"        , emit: gtf
    path "versions.yml" , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${gff.baseName}${task.ext.suffix}" : "${gff.baseName}"
    def args     = task.ext.args ?: ''
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
