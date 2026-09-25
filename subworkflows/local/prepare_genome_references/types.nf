// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
// Superseded or incidental reference files, still published under --save_reference.
record GenomeIntermediates {
    gff:                          Path?
    additional_fasta:             Path?
    gtf_pre_filter:               Path?
    fasta_pre_concat:             Path?
    transcript_fasta_pre_gencode: Path?
    transcript_fasta_rsem_dir:    Path?
}

record GenomeReferences {
    fasta:            Path?
    fai:              Path?
    gtf:              Path?
    gene_bed:         Path?
    transcript_fasta: Path?
    chrom_sizes:      Path?
    rrna_fastas:      List<Path>?
    kraken_db:        Path?
    intermediates:    GenomeIntermediates?
}
