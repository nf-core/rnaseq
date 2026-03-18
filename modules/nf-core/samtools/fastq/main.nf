process SAMTOOLS_FASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(input)
    val interleave

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz"), optional: true, emit: fastq
    tuple val(meta), path("*_interleaved.fastq"), optional: true, emit: interleaved
    tuple val(meta), path("*_singleton.fastq.gz"), optional: true, emit: singleton
    tuple val(meta), path("*_other.fastq.gz"), optional: true, emit: other
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = interleave && !meta.single_end
        ? "> ${prefix}_interleaved.fastq"
        : meta.single_end
            ? "-1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz"
            : "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        fastq \\
        ${args} \\
        --threads ${task.cpus - 1} \\
        -0 ${prefix}_other.fastq.gz \\
        ${input} \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = interleave && !meta.single_end
        ? "touch ${prefix}_interleaved.fastq"
        : meta.single_end
            ? "echo | gzip > ${prefix}_1.fastq.gz && echo | gzip > ${prefix}_singleton.fastq.gz"
            : "echo | gzip > ${prefix}_1.fastq.gz && echo | gzip > ${prefix}_2.fastq.gz && echo | gzip > ${prefix}_singleton.fastq.gz"
    """
    ${output}
    echo | gzip > ${prefix}_other.fastq.gz
    """
}
