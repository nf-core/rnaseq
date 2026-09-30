nextflow.enable.types = true

include { ReadsInput } from '../../types'

process SYLPH_PROFILE {
    tag "${sample.meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    sample: ReadsInput
    database: List<Path>

    output:
    record(id: sample.id, meta: sample.meta, profile_out: file('*.tsv'))

    topic:
    tuple(task.process, 'sylph', eval('sylph -V | sed "s/sylph //g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def input = sample.meta.single_end ? "-r ${sample.reads[0]}" : "-1 ${sample.reads[0]} -2 ${sample.reads[1]}"
    """
    sylph profile \\
        -t ${task.cpus} \\
        ${args} \\
        ${database.join(' ')} \\
        ${input} \\
        -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.tsv
    """
}
