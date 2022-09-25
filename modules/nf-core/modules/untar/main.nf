process UNTAR {
    tag "$archive"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$untar"), emit: untar
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    untar     = archive.toString() - '.tar.gz'

    """
    mkdir output

    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in output
    if [[ \$(tar -tzf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C output --strip-components 1 \\
            -xzvf \\
            $args \\
            $archive \\
            $args2
    else
        tar \\
            -C output \\
            -xzvf \\
            $args \\
            $archive \\
            $args2
    fi

    mv output ${untar}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    untar     = archive.toString() - '.tar.gz'
    """
    touch $untar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
