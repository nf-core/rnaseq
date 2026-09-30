// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

// rsem_merge emit: record(id, rsem_merge: RsemMerge)
record RsemMerge {
    counts_gene:       Path
    tpm_gene:          Path
    counts_transcript: Path
    tpm_transcript:    Path
    genes_long:        Path
    isoforms_long:     Path
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
}
