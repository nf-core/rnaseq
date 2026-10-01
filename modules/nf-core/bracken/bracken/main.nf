nextflow.enable.types = true

include { BrackenResult } from '../../types'

record BrackenInput {
    id:     String
    meta:   Map
    report: Path
    args:   String?
}

process BRACKEN_BRACKEN {
    tag "${sample.meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data':
        'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01' }"

    input:
    sample: BrackenInput
    database: Path

    output:
    record(
        id:        sample.id,
        meta:      sample.meta,
        abundance: file(bracken_report),
        report:    file(bracken_kraken_style_report)
    ) as BrackenResult

    topic:
    tuple(task.process, 'bracken', eval('bracken -v | cut -f2 -d"v"')) >> 'versions'

    script:
    def args = [sample.args, task.ext.args].findAll { a -> a }.join(' ')
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${sample.report}' \\
        -o '${bracken_report}' \\
        -w '${bracken_kraken_style_report}'
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    """
    touch ${prefix}.tsv
    touch ${bracken_kraken_style_report}
    """
}
