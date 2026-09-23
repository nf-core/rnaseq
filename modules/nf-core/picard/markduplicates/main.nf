process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b474c6f12c0502377b95062b75ef4b7778864b1353c7fc1d5a3b1f3b3017fd2e/data'
        : 'community.wave.seqera.io/library/picard:3.5.0--842d4c70c98af9b4'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('picard'), eval("picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = reads.getExtension()
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    if ("${reads}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        ${args} \\
        --INPUT ${reads} \\
        --OUTPUT ${prefix}.${suffix} \\
        ${reference} \\
        --METRICS_FILE ${prefix}.metrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = reads.getExtension()
    if ("${reads}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.${suffix}.bai
    touch ${prefix}.metrics.txt
    """
}
