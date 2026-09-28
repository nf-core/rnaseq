// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record StringtieAssembly {
    id:             String
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
}

record StringtieMerged {
    id:         String
    merged_gtf: Path
}
