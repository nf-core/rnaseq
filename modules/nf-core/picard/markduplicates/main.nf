nextflow.preview.types = true

record PicardMarkDupResult {
    meta:    Map
    bam:     Path?
    bai:     Path?
    cram:    Path?
    metrics: Path
}

process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    (meta: Map, reads: Path): Record
    (meta2: Map, fasta: Path): Record
    (meta3: Map, fai: Path): Record

    output:
    record(
        meta:    meta,
        bam:     file("*.bam",          optional: true),
        bai:     file("*.bai",          optional: true),
        cram:    file("*.cram",         optional: true),
        metrics: file("*.metrics.txt")
    )
    tuple val("${task.process}"), val('picard'), eval('echo $(picard MarkDuplicates --version 2>&1) | grep -o "Version:.*" | cut -f2- -d:'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $reads \\
        --OUTPUT ${prefix}.${suffix} \\
        $reference \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.MarkDuplicates.metrics.txt
    """
}
