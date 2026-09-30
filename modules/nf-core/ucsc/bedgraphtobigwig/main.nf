nextflow.enable.types = true

include { BedgraphInput } from '../../types'

record UcscBedgraphtobigwigResult {
    id:     String
    meta:   Map
    bigwig: Path
}

process UCSC_BEDGRAPHTOBIGWIG {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:482--hdc0a859_0' :
        'quay.io/biocontainers/ucsc-bedgraphtobigwig:482--hdc0a859_0' }"

    input:
    sample: BedgraphInput
    sizes: Path

    output:
    record(id: sample.id, meta: sample.meta, bigwig: file("*.bigWig"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'ucsc', '482') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    bedGraphToBigWig \\
        $args \\
        ${sample.bedgraph} \\
        $sizes \\
        ${prefix}.bigWig
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.bigWig
    """
}
