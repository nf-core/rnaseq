process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda "bioconda::gffread=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    path gff

    output:
    path "*.gtf"        , emit: gtf
    tuple val("${task.process}"), val('gffread'), cmd("gffread --version 2>&1"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${gff.baseName}"
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.gtf
    """
}
