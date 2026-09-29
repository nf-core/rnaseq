nextflow.enable.types = true

process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data':
        'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01' }"

    input:
    tuple(meta: Map, kraken_report: Path)
    database: Path

    output:
    record(
        meta:      meta,
        abundance: file("${task.ext.prefix ?: meta.id}.tsv"),
        report:    file("${task.ext.prefix ?: meta.id}.kraken2.report_bracken.txt")
    )

    topic:
    tuple(task.process, 'bracken', eval('bracken -v | cut -f2 -d"v"')) >> 'versions'

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${kraken_report}' \\
        -o '${prefix}.tsv' \\
        -w '${prefix}.kraken2.report_bracken.txt'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.kraken2.report_bracken.txt
    """
}
