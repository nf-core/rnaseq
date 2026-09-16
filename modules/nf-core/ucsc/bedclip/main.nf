process UCSC_BEDCLIP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedclip:377--h0b8a92a_2' :
        'biocontainers/ucsc-bedclip:377--h0b8a92a_2' }"

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val("${task.process}"), val('ucsc'), eval("echo $VERSION"), topic: versions, emit: versions_ucsc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bedClip \\
        $args \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bedGraph
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bedGraph
    """
}
