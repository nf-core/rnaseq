process HISAT2_EXTRACTSPLICESITES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f546b522b78597047a54609eefcf35a014de082df0e500bc2adf64d2733669/data'
        : 'community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.splice_sites.txt"), emit: txt
    tuple val("${task.process}"), val('hisat2'), eval('hisat2 --version | grep -o "version [^ ]*" | cut -d " " -f 2'), topic: versions, emit: versions_hisat2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py \\
        ${args} \\
        ${gtf} \\
        > ${gtf.baseName}.splice_sites.txt
    """

    stub:
    """
    touch ${gtf.baseName}.splice_sites.txt
    """
}
