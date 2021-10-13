process UNTAR {
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
    path "$untar"      , emit: untar
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    untar    = archive.toString() - '.tar.gz'
    """
    tar \\
        -xzvf \\
        $args \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
