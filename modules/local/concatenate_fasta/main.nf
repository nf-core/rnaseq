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

    input:
    sample: ConcatenateFastaInput

    output:
    record(id: sample.id, meta: sample.meta, fasta: file('rrna_combined_dna.fasta')) as ConcatenateFastaResult

    exec:
    def combined = task.workDir.resolve('rrna_combined_dna.fasta')
    sample.fastas.each { f ->
        combined << f.text
        combined << '\n'
    }
}
