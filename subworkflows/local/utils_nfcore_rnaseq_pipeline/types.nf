// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).

// One independently-produced, optionally-user-supplied publish artifact. Used to stream
// genome reference/index/intermediate files as records instead of fusing them into one
// wide row via a fake join - each producer mixes in its own record, or none at all.
record GenomeArtifact {
    kind: String
    file: Path
}
