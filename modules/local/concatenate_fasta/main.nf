nextflow.enable.types = true

record ConcatenateFastaInput {
    id:     String
    meta:   Map
    fastas: List<Path>
}

record ConcatenateFastaResult {
    id:    String
    meta:  Map
    fasta: Path
}

process CONCATENATE_FASTA {
    tag "${sample.meta.id}"
    label 'process_single'

    input:
    sample: ConcatenateFastaInput

    output:
    record(id: sample.id, meta: sample.meta, fasta: file('rrna_combined_dna.fasta')) as ConcatenateFastaResult

    script:
    def fastas = sample.fastas.collect { f -> f.name }.join(' ')
    """
    for f in ${fastas}; do
        cat "\$f"
        echo
    done > rrna_combined_dna.fasta
    """

    stub:
    """
    touch rrna_combined_dna.fasta
    """
}
