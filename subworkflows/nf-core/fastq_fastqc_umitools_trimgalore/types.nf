record TrimgaloreFastqc {
    raw_html: List<Path>
    raw_zip:  List<Path>
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

record FastqFastqcUmitoolsTrimgalore {
    id:                String
    fastqc:            TrimgaloreFastqc?
    umi:               UmitoolsExtractFiles?
    trim:              TrimgaloreTrim?
    num_trimmed_reads: Float?
}
