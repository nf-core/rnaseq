// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record FastqRemoveRrna {
    id:               String
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
}
