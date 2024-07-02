process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff)
    path fasta

    output:
    tuple val(meta), path("*.gtf")  , emit: gtf             , optional: true
    tuple val(meta), path("*.gff3") , emit: gffread_gff     , optional: true
    tuple val(meta), path("*.fasta"), emit: gffread_fasta   , optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { args.contains(it) } ) ? 'fasta' : 'gff3' )
    def fasta_arg   = fasta                     ? "-g $fasta" : ''
    def output_name = "${prefix}.${extension}"
    def output      = extension == "fasta"      ? "$output_name" : "-o $output_name"
    def args_sorted = args.replaceAll(/(.*)(-[wxy])(.*)/) { all, pre, param, post -> "$pre $post $param" }.trim()
    // args_sorted  = Move '-w', '-x', and '-y' to the end of the args string as gffread expects the file name after these parameters
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    gffread \\
        $gff \\
        $fasta_arg \\
        $args_sorted \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { args.contains(it) } ) ? 'fasta' : 'gff3' )
    def output_name = "${prefix}.${extension}"
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch $output_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
