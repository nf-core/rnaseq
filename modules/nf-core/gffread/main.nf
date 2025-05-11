process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    gff     : Path

    output:
    gtf         : Path? = file("*.gtf")
    gffread_gff : Path? = file("*.gff3")

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${gff.baseName}"
    def extension   = args.contains("-T") ? 'gtf' : 'gffread.gff3'
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.${extension}
    """
}
