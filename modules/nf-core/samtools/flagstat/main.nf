process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        ${bam} \\
        > ${prefix}.flagstat
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_FLAGSTAT > ${prefix}.flagstat
    1000000 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    900000 + 0 mapped (90.00% : N/A)
    1000000 + 0 paired in sequencing
    500000 + 0 read1
    500000 + 0 read2
    800000 + 0 properly paired (80.00% : N/A)
    850000 + 0 with mate mapped to a different chr
    50000 + 0 with mate mapped to a different chr (mapQ>=5)
    END_FLAGSTAT
    """
}
