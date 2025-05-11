process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    meta    : Map
    bam     : Path
    bai     : Path

    output:
    file("*.idxstats")

    topic:
    file('versions.yml') >> 'versions'
    file("*.idxstats") >> 'logs'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        idxstats \\
        --threads ${task.cpus-1} \\
        $bam \\
        > ${prefix}.idxstats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.idxstats
    """
}
