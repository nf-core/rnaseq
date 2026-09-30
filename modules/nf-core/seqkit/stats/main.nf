nextflow.enable.types = true

process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    record(id: String, meta: Map, reads: List<Path>)

    output:
    record(id: id, meta: meta, stats: file("*.tsv"))

    topic:
    tuple(task.process, 'seqkit', eval("seqkit version | sed 's/^.*v//'")) >> 'versions'

    script:
    def args = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads.join(' ')} > '${prefix}.tsv'
    """

    stub:
    def args = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.tsv
    """
}
