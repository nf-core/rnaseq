// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

record TrimgaloreTrim {
    html:     List<Path>?
    zip:      List<Path>?
    log:      List<Path>?
    json:     List<Path>?
    unpaired: List<Path>?
}

record UmitoolsExtractFiles {
    log:   Path
    reads: List<Path>
}

// `reads` and `meta` are the reads handed to the next stage: UMI-extracted and
// with R1 or R2 discarded when requested, then trimmed. `reads` is null for
// samples below min_trimmed_reads. TrimgaloreTrim.html and zip are FastQC
// reports on the trimmed reads.
record FastqFastqcUmitoolsTrimgalore {
    id:                String
    meta:              Map
    reads:             List<Path>?
    fastqc_raw_html:   List<Path>?
    fastqc_raw_zip:    List<Path>?
    umi:               UmitoolsExtractFiles?
    trim:              TrimgaloreTrim?
    num_trimmed_reads: Float?
}
