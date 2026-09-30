nextflow.enable.types = true

include { BamInput } from '../../types'

process SAMTOOLS_FASTQ {
    tag "${sample.meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: BamInput
    interleave: Boolean

    output:
    record(
        id:          sample.id,
        meta:        sample.meta,
        reads:       files('*_{1,2}.fastq.gz', optional: true).toSorted { f -> f.name },
        interleaved: file('*_interleaved.fastq', optional: true),
        singleton:   file('*_singleton.fastq.gz', optional: true),
        other:       file('*_other.fastq.gz', optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def output = interleave && !sample.meta.single_end
        ? "> ${prefix}_interleaved.fastq"
        : sample.meta.single_end
            ? "-1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz"
            : "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        fastq \\
        ${args} \\
        --threads ${task.cpus - 1} \\
        -0 ${prefix}_other.fastq.gz \\
        ${sample.bam} \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def output_command = interleave && !sample.meta.single_end
        ? "touch ${prefix}_interleaved.fastq"
        : sample.meta.single_end
            ? "echo | bgzip -c > ${prefix}_1.fastq.gz && echo | bgzip -c > ${prefix}_singleton.fastq.gz"
            : "echo | bgzip -c > ${prefix}_1.fastq.gz && echo | bgzip -c > ${prefix}_2.fastq.gz && echo | bgzip -c > ${prefix}_singleton.fastq.gz"
    def other_command = "echo | bgzip -c > ${prefix}_other.fastq.gz"
    """
    ${output_command}
    ${other_command}
    """
}
