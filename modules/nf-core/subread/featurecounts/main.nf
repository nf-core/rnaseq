nextflow.enable.types = true

process SUBREAD_FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0'
        : 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'}"

    input:
    tuple(meta: Map, bams: List<Path>, annotation: Path)

    output:
    record(
        meta:    meta,
        counts:  file("*featureCounts.tsv"),
        summary: file("*featureCounts.tsv.summary")
    )

    topic:
    tuple(task.process, 'subread', eval("featureCounts -v 2>&1 | sed 's/featureCounts v//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    }
    else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    """
    featureCounts \\
        ${args} \\
        ${paired_end} \\
        -T ${task.cpus} \\
        -a ${annotation} \\
        -s ${strandedness} \\
        -o ${prefix}.featureCounts.tsv \\
        ${bams.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    """
}
