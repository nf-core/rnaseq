nextflow.preview.types = true

record SamtoolsIndexResult {
    meta: Map
    bai:  Path?
    csi:  Path?
    crai: Path?
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    (meta: Map, input: Path): Record

    output:
    record(
        meta: meta,
        bai:  file("*.bai",  optional: true),
        csi:  file("*.csi",  optional: true),
        crai: file("*.crai", optional: true)
    )
    tuple val("${task.process}"), val('samtools'), eval('echo $(samtools --version 2>&1) | sed "s/^.*samtools //; s/Using.*$//"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        $args \\
        $input
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram' ?
                    "crai" : args.contains("-c") ?  "csi" : "bai"
    """
    touch ${input}.${extension}
    """
}
