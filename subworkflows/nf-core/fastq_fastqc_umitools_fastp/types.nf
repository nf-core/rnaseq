// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record FastpReads {
    id:            String
    meta:          Map
    reads:         List<Path>
    adapter_fasta: Path?
}

record FastpTrim {
    html:         Path
    json:         Path
    log:          Path
    reads_fail:   List<Path>?
    reads_merged: Path?
}

record UmitoolsExtractFiles {
    log:   Path
    reads: List<Path>
}

// `reads` and `meta` are the reads handed to the next stage: UMI-extracted and
// with R1 or R2 discarded when requested, then trimmed. `reads` is null for
// samples below min_trimmed_reads or without any trimmed reads.
record FastqFastqcUmitoolsFastp {
    id:                String
    meta:              Map
    reads:             List<Path>?
    fastqc_raw_html:   List<Path>?
    fastqc_raw_zip:    List<Path>?
    fastqc_trim_html:  List<Path>?
    fastqc_trim_zip:   List<Path>?
    umi:               UmitoolsExtractFiles?
    trim:              FastpTrim?
    adapter_seq:       String?
    num_trimmed_reads: Long?
}
