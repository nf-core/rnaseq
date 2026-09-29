nextflow.enable.types = true

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.2.0--hdcf5f25_6':
        'quay.io/biocontainers/preseq:3.2.0--hdcf5f25_6' }"

    input:
    tuple(meta: Map, bam: Path)

    output:
    record(
        meta:      meta,
        lc_extrap: file("*.lc_extrap.txt"),
        log:       file("*.log")
    )

    topic:
    tuple(task.process, 'preseq', eval("preseq 2>&1 | sed -n 's/Version: //p'")) >> 'versions'

    script:
    def args = (task.ext.args ?: '') + (task.attempt > 1 ? ' -defects' : '')  // Disable testing for defects
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        ${args} \\
        ${paired_end} \\
        -output ${prefix}.lc_extrap.txt \\
        ${bam}
    cp .command.err ${prefix}.command.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lc_extrap.txt
    touch ${prefix}.command.log
    """
}
