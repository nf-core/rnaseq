process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    path gff

    output:
    path "*.gtf"        , emit: gtf         , optional: true
    path "*.gff3"       , emit: gffread_gff , optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${gff.baseName}"
    def extension   = args.contains("-T") ? 'gtf' : 'gffread.gff3'
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.${extension}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
