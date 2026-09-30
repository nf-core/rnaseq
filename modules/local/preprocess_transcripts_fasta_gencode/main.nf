nextflow.enable.types = true

include { FastaInput } from '../../nf-core/types'

process PREPROCESS_TRANSCRIPTS_FASTA_GENCODE {
    tag "${sample.fasta}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    sample: FastaInput

    output:
    record(id: sample.id, meta: sample.meta, fasta: file('*.fa'))

    topic:
    tuple(task.process, 'sed', eval("sed --version 2>&1 | sed '1!d;s/^.*) //'")) >> 'versions'

    script:
    def gzipped = sample.fasta.name.endsWith('.gz')
    def outfile = gzipped ? file(sample.fasta.baseName).baseName : sample.fasta.baseName
    def command = gzipped ? 'zcat' : 'cat'
    """
    $command ${sample.fasta} | cut -d "|" -f1 > ${outfile}.fixed.fa
    """

    stub:
    def gzipped = sample.fasta.name.endsWith('.gz')
    def outfile = gzipped ? file(sample.fasta.baseName).baseName : sample.fasta.baseName
    """
    touch ${outfile}.fixed.fa
    """
}
