nextflow.enable.types = true


record SamtoolsIndexInput {
    id:   String
    meta: Map
    bam:  Path
    args: String?
}

record SamtoolsIndexResult {
    id:   String
    meta: Map
    bai:  Path
}

process SAMTOOLS_INDEX {
    tag "${sample.meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: SamtoolsIndexInput

    output:
    record(id: sample.id, meta: sample.meta, bai: file('*.{bai,csi,crai}')) as SamtoolsIndexResult

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: sample.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${sample.bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = sample.bam.getExtension() == 'cram'
        ? "crai"
        : args.contains("-c") ? "csi" : "bai"
    """
    touch ${sample.bam}.${extension}
    """
}
