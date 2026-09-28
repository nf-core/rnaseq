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
