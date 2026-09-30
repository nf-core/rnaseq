nextflow.enable.types = true

process SYLPH_PROFILE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    record(id: String, meta: Map, reads: List<Path>)
    database: List<Path>

    output:
    record(id: id, meta: meta, profile_out: file('*.tsv'))

    topic:
    tuple(task.process, 'sylph', eval('sylph -V | sed "s/sylph //g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph profile \\
        -t ${task.cpus} \\
        ${args} \\
        ${database.join(' ')} \\
        ${input} \\
        -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
