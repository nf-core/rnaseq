// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record BigwigFiles {
    id:       String
    bigwig:   Path
    bedgraph: Path
}
