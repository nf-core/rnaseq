record RseqcInnerDistance {
    distance: Path
    freq:     Path?
    mean:     Path?
    pdf:      Path?
    rscript:  Path?
}

record RseqcJunctionAnnotation {
    bed:          Path?
    interact_bed: Path?
    xls:          Path
    pdf:          Path?
    events_pdf:   Path?
    rscript:      Path
    log:          Path
}

record RseqcJunctionSaturation {
    pdf:     Path
    rscript: Path
}

record RseqcReadDuplication {
    seq_xls: Path
    pos_xls: Path
    pdf:     Path
    rscript: Path
}

record RseqcTin {
    txt: Path
    xls: Path
}

record Rseqc {
    id:                 String
    bamstat:            Path?
    inferexperiment:    Path?
    innerdistance:      RseqcInnerDistance?
    junctionannotation: RseqcJunctionAnnotation?
    junctionsaturation: RseqcJunctionSaturation?
    readdistribution:   Path?
    readduplication:    RseqcReadDuplication?
    tin:                RseqcTin?
}
