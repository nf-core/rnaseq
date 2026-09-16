process EAUTILS_GTF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed") , emit: bed
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtf2bed.pl'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V, 1)')
    END_VERSIONS
    """
}
