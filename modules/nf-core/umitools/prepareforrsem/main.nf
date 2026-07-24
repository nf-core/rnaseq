process UMITOOLS_PREPAREFORRSEM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
        'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log
    tuple val("${task.process}"), val('umitools'), eval("umi_tools --version | sed 's/UMI-tools version: //'"), emit: versions_umitools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    umi_tools prepare-for-rsem \\
        --stdin=$bam \\
        --stdout=${prefix}.bam \\
        --log=${prefix}.prepare_for_rsem.log \\
        $args
    """

    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.log
    """
}
