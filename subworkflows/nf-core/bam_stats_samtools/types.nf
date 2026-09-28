// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record SamtoolsStats {
    id:       String
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}
