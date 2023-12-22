process PREPROCESS_TRANSCRIPTS_FASTA_GENCODE {
    tag "$fasta"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path fasta

    output:
    path "*.fa"        , emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gzipped = fasta.toString().endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    def command = gzipped ? 'zcat' : 'cat'
    """
    $command $fasta | cut -d "|" -f1 > ${outfile}.fixed.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
