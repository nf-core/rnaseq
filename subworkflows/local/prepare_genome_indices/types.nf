// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record GenomeIndices {
    star:                  Path?
    rsem:                  Path?
    rsem_transcript_fasta: Path?
    hisat2:                Path?
    hisat2_splicesites:    Path?
    bowtie2:               Path?
    salmon:                Path?
    kallisto:              Path?
    bbsplit:               Path?
    bbsplit_log:           Path?
    sortmerna:             Path?
    bowtie2_rrna:          Path?
}
