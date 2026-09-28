record GenomeReferences {
    fasta:            Path?
    fai:              Path?
    gtf:              Path?
    gene_bed:         Path?
    transcript_fasta: Path?
    chrom_sizes:      Path?
    rrna_fastas:      List<Path>?
    kraken_db:        Path?
}
