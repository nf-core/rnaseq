// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Bedgraph {
    id:       String
    meta:     Map
    bedgraph: Path
}

record BigwigFiles {
    id:       String
    meta:     Map
    bedgraph: Path
    bigwig:   Path
}
