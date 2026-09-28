record FastpFastqc {
    raw_html:  List<Path>
    raw_zip:   List<Path>
    trim_html: List<Path>?
    trim_zip:  List<Path>?
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

record FastqFastqcUmitoolsFastp {
    id:                String
    fastqc:            FastpFastqc?
    umi:               UmitoolsExtractFiles?
    trim:              FastpTrim?
    adapter_seq:       String?
    num_trimmed_reads: Long?
}
