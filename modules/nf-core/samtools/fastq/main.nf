process SAMTOOLS_FASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

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
    def output_command = interleave && !meta.single_end
        ? "touch ${prefix}_interleaved.fastq"
        : meta.single_end
            ? "echo | bgzip -c > ${prefix}_1.fastq.gz && echo | bgzip -c > ${prefix}_singleton.fastq.gz"
            : "echo | bgzip -c > ${prefix}_1.fastq.gz && echo | bgzip -c > ${prefix}_2.fastq.gz && echo | bgzip -c > ${prefix}_singleton.fastq.gz"
    def other_command = "echo | bgzip -c > ${prefix}_other.fastq.gz"
    """
    ${output_command}
    ${other_command}
    """
}
