// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

record SalmonQuantSample {
    id:                String
    meta:              Map
    quant_dir:         Path
    json_info:         Path?
    lib_format_counts: Path?
}

record KallistoQuantSample {
    id:        String
    meta:      Map
    quant_dir: Path
    json_info: Path
    log:       Path
}
