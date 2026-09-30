// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

// `reads` is null for samples with no reads left after rRNA removal. The tool
// logs are set only for the selected tool.
record FastqRemoveRrna {
    id:               String
    meta:             Map
    reads:            List<Path>?
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
}

// Run-level references, each set only when it is built here.
record RrnaReferences {
    sortmerna_index:  Path?
    bowtie2_index:    Path?
    seqkit_prefixed:  List<Path>?
    seqkit_converted: List<Path>?
}
