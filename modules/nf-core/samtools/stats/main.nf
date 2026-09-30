nextflow.enable.types = true

include { BamBaiInput } from '../../types'

record SamtoolsStatsResult {
    id:    String
    meta:  Map
    stats: Path
}

process SAMTOOLS_STATS {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: BamBaiInput
    fasta: Path?
    fai: Path?

    output:
    record(id: sample.id, meta: sample.meta, stats: file('*.stats'))

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${sample.bam} \\
        > ${prefix}.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.stats
    """
}
