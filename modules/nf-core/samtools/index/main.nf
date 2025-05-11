process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    meta    : Map
    input   : Path

    output:
    bai     : Path = file("*.bai") 
    csi     : Path = file("*.csi") 
    crai    : Path = file("*.crai")

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input
    """

    stub:
    """
    touch ${input}.bai
    touch ${input}.crai
    touch ${input}.csi
    """
}
