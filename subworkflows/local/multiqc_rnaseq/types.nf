// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record MultiqcReport {
    id:     String
    meta:   Map
    report: Path
    data:   Path
    plots:  Path?
}
