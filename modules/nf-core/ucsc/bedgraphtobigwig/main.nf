process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:469--h9b8f530_0' :
        'biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0' }"

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    tuple val("${task.process}"), val('ucsc'), eval("echo $VERSION"), topic: versions, emit: versions_ucsc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '469' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bedGraphToBigWig \\
        $args \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bigWig
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '469' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bigWig
    """
}
