include { FastqRemoveRrna } from '../fastq_remove_rrna/types'

// TrimGalore and fastp report the same kinds of file with different
// cardinality, so every multi-file field is a list and tool-specific fields
// are nullable. TrimGalore's own html/zip are FastQC reports on the trimmed
// reads and are held in fastqc.trim_html/trim_zip.
record PreprocessedFastqc {
    raw_html:      List<Path>?
    raw_zip:       List<Path>?
    trim_html:     List<Path>?
    trim_zip:      List<Path>?
    filtered_html: List<Path>?
    filtered_zip:  List<Path>?
}

record PreprocessedTrim {
    html:         List<Path>?
    log:          List<Path>?
    json:         List<Path>?
    unpaired:     List<Path>?
    reads_fail:   List<Path>?
    reads_merged: Path?
}

record PreprocessedUmi {
    log:   Path
    reads: List<Path>
}

record PreprocessedBbsplit {
    stats: Path
}

record PreprocessedLint {
    raw:     Path?
    trimmed: Path?
    bbsplit: Path?
    ribo:    Path?
}

// Samples that fail min_trimmed_reads keep their record with null reads and
// reads_trimmed, and a meta that lacks the inferred strandedness.
record FastqQcTrimFilterSetstrandedness {
    id:                String
    meta:              Map
    reads:             List<Path>?
    reads_cat:         List<Path>
    reads_trimmed:     List<Path>?
    num_trimmed_reads: Long?
    fastqc:            PreprocessedFastqc?
    trim:              PreprocessedTrim?
    umi:               PreprocessedUmi?
    bbsplit:           PreprocessedBbsplit?
    lint:              PreprocessedLint?
    rrna:              FastqRemoveRrna?
}
