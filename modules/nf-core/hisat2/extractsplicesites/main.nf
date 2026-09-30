nextflow.enable.types = true

include { GtfInput } from '../../types'

process HISAT2_EXTRACTSPLICESITES {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f546b522b78597047a54609eefcf35a014de082df0e500bc2adf64d2733669/data'
        : 'community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89'}"

    input:
    sample: GtfInput

    output:
    record(id: sample.id, meta: sample.meta, splicesites: file('*.splice_sites.txt'))

    topic:
    tuple(task.process, 'hisat2', eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2')) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools --version | sed -n '1s/samtools //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py \\
        ${args} \\
        ${sample.gtf} \\
        > ${sample.gtf.baseName}.splice_sites.txt
    """

    stub:
    """
    touch ${sample.gtf.baseName}.splice_sites.txt
    """
}
