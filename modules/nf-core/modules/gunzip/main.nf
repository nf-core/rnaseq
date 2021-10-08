// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

process GUNZIP {
    tag "$archive"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path archive

    output:
    path "$gunzip",       emit: gunzip
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    gunzip   = archive.toString() - '.gz'
    """
    gunzip \\
        -f \\
        $args \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
