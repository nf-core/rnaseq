nextflow.enable.types = true

process UCSC_BEDCLIP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedclip:482--h0b57e2e_0' :
        'quay.io/biocontainers/ucsc-bedclip:482--h0b57e2e_0' }"

    input:
    tuple(meta: Map, bedgraph: Path)
    sizes: Path

    output:
    record(meta: meta, bedgraph: file("*.bedGraph"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'ucsc', '482') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedClip \\
        $args \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bedGraph
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedGraph
    """
}
