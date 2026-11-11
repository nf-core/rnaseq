process PREPROCESS_TRANSCRIPTS_FASTA_GENCODE {
    tag "$fasta"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path fasta

    output:
    path "*.fa", emit: fasta
    tuple val("${task.process}"), val('sed'), eval("sed --version 2>&1 | sed '1!d;s/^.*) //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gzipped = fasta.toString().endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    def command = gzipped ? 'zcat' : 'cat'
    """
    $command $fasta | cut -d "|" -f1 > ${outfile}.fixed.fa
    """

    stub:
    def gzipped = fasta.toString().endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    """
    touch ${outfile}.fixed.fa
    """
}
