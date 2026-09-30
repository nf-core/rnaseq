nextflow.enable.types = true

include { BamBaiInput } from '../../types'

record SamtoolsIdxstatsResult {
    id:       String
    meta:     Map
    idxstats: Path
}

process SAMTOOLS_IDXSTATS {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: BamBaiInput

    output:
    record(id: sample.id, meta: sample.meta, idxstats: file('*.idxstats'))

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"

    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        idxstats \\
        --threads ${task.cpus - 1} \\
        ${sample.bam} \\
        > ${prefix}.idxstats
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"

    """
    touch ${prefix}.idxstats
    """
}
