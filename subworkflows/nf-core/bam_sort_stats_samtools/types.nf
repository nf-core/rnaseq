// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record RawBams {
    id:       String
    meta:     Map
    raw_bams: List<Path>
}
