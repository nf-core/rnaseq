nextflow.preview.types = true

record PrepareForRsemResult {
    meta: Map
    bam:  Path
    log:  Path
}

process UMITOOLS_PREPAREFORRSEM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.6--py311haab0aaa_0' :
        'biocontainers/umi_tools:1.1.6--py311haab0aaa_0' }"

    input:
    (meta: Map, bam: Path, bai: Path): Record

    output:
    record(
        meta: meta,
        bam:  file('*.bam'),
        log:  file('*.log')
    )
    tuple val("${task.process}"), val('umi_tools'), eval('umi_tools --version | sed "/version:/!d; s/.*: //"'), topic: versions

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
