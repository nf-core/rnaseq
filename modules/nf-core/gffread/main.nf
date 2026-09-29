nextflow.enable.types = true

process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple(meta: Map, gff: Path)
    fasta: Path?

    output:
    record(
        meta:          meta,
        gtf:           file('*.gtf', optional: true),
        gffread_gff:   file('*.gff3', optional: true),
        gffread_fasta: file('*.fasta', optional: true),
        bed:           file('*.bed', optional: true)
    )

    topic:
    tuple(task.process, 'gffread', eval('gffread --version 2>&1')) >> 'versions'

    script:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def fasta_arg   = fasta                     ? "-g $fasta" : ''
    def output_name = "${prefix}.${extension}"
    def output      = extension == "fasta"      ? "$output_name" : "-o $output_name"
    def args_sorted = args.replaceAll(/(.*)(-[wxy])(.*)/, '$1 $3 $2').trim()
    // args_sorted  = Move '-w', '-x', and '-y' to the end of the args string as gffread expects the file name after these parameters
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    gffread \\
        $gff \\
        $fasta_arg \\
        $args_sorted \\
        $output
    """

    stub:
    def args        = task.ext.args             ?: ''
    def prefix      = task.ext.prefix           ?: "${meta.id}"
    def extension   = args.contains("--bed")    ? 'bed' : ( args.contains("-T")       ? 'gtf' : ( ( ['-w', '-x', '-y' ].any { flag -> args.contains(flag) } ) ? 'fasta' : 'gff3' ) )
    def output_name = "${prefix}.${extension}"
    if ( "$output_name" in [ "$gff", "$fasta" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch $output_name
    """
}
