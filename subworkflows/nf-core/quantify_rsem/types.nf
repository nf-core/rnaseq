record RsemMerge {
    counts_gene:       Path
    tpm_gene:          Path
    counts_transcript: Path
    tpm_transcript:    Path
    genes_long:        Path
    isoforms_long:     Path
}

// The fields of QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT's QuantMerged plus
// rsem_merge. Records cannot extend one another, so keep the shared fields in
// sync with QuantMerged.
record RsemQuantMerged {
    id:                        String
    meta:                      Map
    tpm_gene:                  Path
    counts_gene:               Path
    lengths_gene:              Path
    counts_gene_length_scaled: Path
    counts_gene_scaled:        Path
    tpm_transcript:            Path
    counts_transcript:         Path
    lengths_transcript:        Path
    tx2gene:                   Path
    tx2gene_augmented:         Path
    merged_gene_rds:           Path?
    merged_transcript_rds:     Path?
    rsem_merge:                RsemMerge?
}

record RsemQuantSample {
    id:                String
    meta:              Map
    counts_gene:       Path
    counts_transcript: Path
    stat:              Path
    log:               Path?
    bam_star:          Path?
    bam_genome:        Path?
    bam_transcript:    Path?
    quant_merged:      RsemQuantMerged
}
