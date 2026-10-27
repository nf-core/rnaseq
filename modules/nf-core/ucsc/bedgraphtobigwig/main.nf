process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:482--hdc0a859_0' :
        'quay.io/biocontainers/ucsc-bedgraphtobigwig:482--hdc0a859_0' }"

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

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
