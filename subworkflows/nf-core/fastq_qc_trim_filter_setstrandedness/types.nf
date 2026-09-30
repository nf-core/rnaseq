// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

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
    stats:              Path
    primary_reads:      List<Path>?
    other_genome_reads: List<Path>?
}

record PreprocessedLint {
    raw:     Path?
    trimmed: Path?
    bbsplit: Path?
    ribo:    Path?
}

// The tool logs are set only for the selected tool.
record PreprocessedRrna {
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
}

// Run-level references, each set only when it is built here.
record RrnaReferences {
    sortmerna_index:  Path?
    bowtie2_index:    Path?
    seqkit_prefixed:  List<Path>?
    seqkit_converted: List<Path>?
}

// Samples that fail min_trimmed_reads keep their record with null reads and
// reads_trimmed, and a meta that lacks the inferred strandedness.
record FastqQcTrimFilterSetstrandedness {
    id:                String
    meta:              Map
    reads:             List<Path>?
    reads_cat:         List<Path>
    reads_trimmed:     List<Path>?
    num_trimmed_reads: Long?  // Float for TrimGalore, Long for fastp
    fastqc:            PreprocessedFastqc?
    trim:              PreprocessedTrim?
    umi:               PreprocessedUmi?
    bbsplit:           PreprocessedBbsplit?
    lint:              PreprocessedLint?
    rrna:              PreprocessedRrna?
}
