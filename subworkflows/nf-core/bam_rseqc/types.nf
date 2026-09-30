// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record RseqcInnerDistance {
    id:       String
    meta:     Map
    distance: Path
    freq:     Path?
    mean:     Path?
    pdf:      Path?
    rscript:  Path?
}

record RseqcJunctionAnnotation {
    id:           String
    meta:         Map
    bed:          Path?
    interact_bed: Path?
    xls:          Path
    pdf:          Path?
    events_pdf:   Path?
    rscript:      Path
    log:          Path
}

record RseqcJunctionSaturation {
    id:      String
    meta:    Map
    pdf:     Path
    rscript: Path
}

record RseqcReadDuplication {
    id:      String
    meta:    Map
    seq_xls: Path
    pos_xls: Path
    pdf:     Path
    rscript: Path
}

record RseqcTin {
    id:   String
    meta: Map
    txt:  Path
    xls:  Path
}

record Rseqc {
    id:                 String
    meta:               Map
    bamstat:            Path?
    inferexperiment:    Path?
    innerdistance:      RseqcInnerDistance?
    junctionannotation: RseqcJunctionAnnotation?
    junctionsaturation: RseqcJunctionSaturation?
    readdistribution:   Path?
    readduplication:    RseqcReadDuplication?
    tin:                RseqcTin?
}
