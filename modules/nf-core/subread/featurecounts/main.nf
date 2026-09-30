nextflow.enable.types = true

include { BamInput } from '../../types'

record SubreadFeaturecountsResult {
    id:      String
    meta:    Map
    counts:  Path
    summary: Path
}

process SUBREAD_FEATURECOUNTS {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0'
        : 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'}"

    input:
    sample: BamInput
    annotation: Path

    output:
    record(
        id:      sample.id,
        meta:    sample.meta,
        counts:  file("*featureCounts.tsv"),
        summary: file("*featureCounts.tsv.summary")
    ) as SubreadFeaturecountsResult

    topic:
    tuple(task.process, 'subread', eval("featureCounts -v 2>&1 | sed 's/featureCounts v//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def paired_end = sample.meta.single_end ? '' : '-p'

    def strandedness = 0
    if (sample.meta.strandedness == 'forward') {
        strandedness = 1
    }
    else if (sample.meta.strandedness == 'reverse') {
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
        ${sample.bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    """
}
