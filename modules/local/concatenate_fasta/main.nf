nextflow.enable.types = true

process CONCATENATE_FASTA {
    tag "${meta.id}"

    input:
    record(id: String, meta: Map, fastas: List<Path>)

    output:
    record(id: id, meta: meta, fasta: file('rrna_combined_dna.fasta'))

    exec:
    def combined = task.workDir.resolve('rrna_combined_dna.fasta')
    fastas.each { f ->
        combined << f.text
        combined << '\n'
    }
}
