// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record BamBai {
    id:   String
    meta: Map
    bam:  Path
    bai:  Path
}

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record SamtoolsStats {
    id:       String
    meta:     Map
    samtools: SamtoolsStatsFiles
}
