nextflow.enable.types = true

process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:482--hdc0a859_0' :
        'quay.io/biocontainers/ucsc-bedgraphtobigwig:482--hdc0a859_0' }"

    input:
    tuple(meta: Map, bedgraph: Path)
    sizes: Path

    output:
    record(meta: meta, bigwig: file("*.bigWig"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'ucsc', '482') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedGraphToBigWig \\
        $args \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bigWig
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bigWig
    """
}
