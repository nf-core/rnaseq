// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

// `reads` is null for samples with no reads left after rRNA removal. The tool
// logs are set only for the selected tool. The last four fields are run-level
// and repeated on every sample: each is set only when the reference is built
// here, and empty otherwise.
record FastqRemoveRrna {
    id:               String
    meta:             Map
    reads:            List<Path>?
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
    sortmerna_index:  Path?
    bowtie2_index:    Path?
    seqkit_prefixed:  List<Path>?
    seqkit_converted: List<Path>?
}
