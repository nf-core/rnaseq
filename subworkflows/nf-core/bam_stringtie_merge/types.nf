// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record StringtieInput {
    id:     String
    meta:   Map
    bam:    Path?
    lrbam:  Path?
}

record StringtieAssembly {
    id:             String
    meta:           Map
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
}

// One record for the whole run: the merged GTF and the per-sample assemblies that fed it.
record StringtieMerged {
    id:         String
    meta:       Map
    gtf:        List<Path>
    assemblies: List<StringtieAssembly>
    merged_gtf: Path
}
