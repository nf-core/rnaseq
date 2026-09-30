nextflow.enable.types = true

process PREPROCESS_TRANSCRIPTS_FASTA_GENCODE {
    tag "$fasta"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    record(id: String, meta: Map, fasta: Path)

    output:
    record(id: id, meta: meta, fasta: file('*.fa'))

    topic:
    tuple(task.process, 'sed', eval("sed --version 2>&1 | sed '1!d;s/^.*) //'")) >> 'versions'

    script:
    def gzipped = fasta.name.endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    def command = gzipped ? 'zcat' : 'cat'
    """
    $command $fasta | cut -d "|" -f1 > ${outfile}.fixed.fa
    """

    stub:
    def gzipped = fasta.name.endsWith('.gz')
    def outfile = gzipped ? file(fasta.baseName).baseName : fasta.baseName
    """
    touch ${outfile}.fixed.fa
    """
}
