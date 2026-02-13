nextflow.preview.types = true

record PreseqResult {
    meta:      Map
    lc_extrap: Path
    log:       Path
}

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.2.0--hdcf5f25_6':
        'biocontainers/preseq:3.2.0--hdcf5f25_6' }"

    input:
    (meta: Map, bam: Path): Record

    output:
    record(
        meta:      meta,
        lc_extrap: file("*.lc_extrap.txt"),
        log:       file("*.log")
    )
    tuple val("${task.process}"), val('preseq'), eval('echo $(preseq 2>&1) | sed "s/^.*Version: //; s/Usage.*$//"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args = task.attempt > 1 ? args.join(' -defects') : args  // Disable testing for defects
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $args \\
        $paired_end \\
        -output ${prefix}.lc_extrap.txt \\
        $bam
    cp .command.err ${prefix}.command.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lc_extrap.txt
    touch ${prefix}.command.log
    """
}
