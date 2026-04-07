process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram'
        ? "crai"
        : args.contains("-c") ? "csi" : "bai"
    """
    touch ${input}.${extension}
    """
}
