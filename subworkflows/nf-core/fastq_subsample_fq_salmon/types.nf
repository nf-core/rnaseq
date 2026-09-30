// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

// reads are the sub-sampled FASTQ files. index_built is set only when the Salmon index was built here.
record SalmonSubsampled {
    id:                String
    meta:              Map
    reads:             List<Path>
    quant_dir:         Path
    json_info:         Path?
    lib_format_counts: Path?
    index_built:       Path?
}
