process SYLPH_PROFILE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path('*.tsv'), emit: profile_out
    tuple val("${task.process}"), val('sylph'), eval('sylph -V | awk "{print \$2}"'), topic: versions, emit: versions_sylph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph profile \\
        -t ${task.cpus} \\
        ${args} \\
        ${database}\\
        ${input} \\
        -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
