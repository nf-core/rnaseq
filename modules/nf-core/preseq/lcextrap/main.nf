process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_single'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.1.2--h445547b_2':
        'biocontainers/preseq:3.1.2--h445547b_2' }"

    input:
    meta    : Map
    bam     : Path

    output:
    lc_extrap   : Path = file("*.lc_extrap.txt")
    log         : Path = file("*.command.log")

    topic:
    file('versions.yml') >> 'versions'
    file('*.lc_extrap.txt') >> 'logs'

    script:
    def args = task.ext.args ?: ''
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
}
