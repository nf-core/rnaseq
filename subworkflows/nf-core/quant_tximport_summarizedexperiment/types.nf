// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record QuantMerged {
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
}
