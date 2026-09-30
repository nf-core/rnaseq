nextflow.enable.types = true

include { BamBaiInput } from '../../types'

record SamtoolsFlagstatResult {
    id:       String
    meta:     Map
    flagstat: Path
}

process SAMTOOLS_FLAGSTAT {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: BamBaiInput

    output:
    record(id: sample.id, meta: sample.meta, flagstat: file('*.flagstat'))

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        ${sample.bam} \\
        > ${prefix}.flagstat
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
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
