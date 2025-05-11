process GUNZIP {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    archive : Path

    output:
    file(gunzip)

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    # Not calling gunzip itself because it creates files
    # with the original group ownership rather than the
    # default one for that user / the work directory
    gzip \\
        -cd \\
        $args \\
        $archive \\
        > $gunzip
    """

    stub:
    gunzip = archive.toString() - '.gz'
    """
    touch $gunzip
    """
}
