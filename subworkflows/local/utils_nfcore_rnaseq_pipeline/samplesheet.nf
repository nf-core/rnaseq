//
// Samplesheet handling. Kept untyped because `nextflow lint` rejects calls to plugin functions such
// as `validate` in typed scripts (nextflow-io/nextflow#7720).
//

include { validate                  } from 'plugin/nf-schema'
include { checkSamplesAfterGrouping } from './main'

//
// Samplesheet columns in schema order. `sample` is stored in the meta map as `id`.
//
def samplesheetColumns() {
    return [ 'sample', 'fastq_1', 'fastq_2', 'strandedness', 'seq_platform', 'seq_center', 'genome_bam', 'transcriptome_bam', 'percent_mapped' ]
}

//
// A samplesheet row as the plain map that JSON-schema validation expects: unset columns are left out
// and files are given as strings, because the validator cannot serialise Path objects.
//
def rowToSchemaMap(row) {
    return samplesheetColumns()
        .findAll { column -> row[column] != null }
        .collectEntries { column ->
            def value = row[column]
            [ column, value instanceof Path ? value.toUriString() : value ]
        }
}

//
// Validate samplesheet rows against the samplesheet schema, with the same checks that are applied
// to a samplesheet file.
//
def validateSamplesheetRows(rows, schema) {
    if( !rows ) {
        error("The input samplesheet contains no samples")
    }
    validate(rows.collect { row -> rowToSchemaMap(row) }, schema)
}

//
// Samplesheet rows as CSV text. Only the columns that are set in at least one row are written, so
// the text has the columns of the samplesheet the rows came from.
//
def samplesheetRowsToCsv(rows) {
    def columns = samplesheetColumns().findAll { column -> rows.any { row -> row[column] != null } }
    def lines = rows.collect { row ->
        def map = rowToSchemaMap(row)
        columns.collect { column -> map[column] != null ? map[column].toString() : '' }.join(',')
    }
    return ( [ columns.join(',') ] + lines.sort(false) ).join('\n') + '\n'
}

//
// Validate the samplesheet rows and merge the rows of a sample that has several sequencing runs
// into one plain map per sample. Rows are ordered by file path first: the channel they come from
// has no fixed order, and the order of the runs decides the order in which they are merged.
// Keys: id, meta, reads (all FASTQ files), runs (one list of FASTQ files per run), bam,
// transcriptome_bam, percent_mapped and prealigned. A sample is prealigned when alignment is
// skipped and the samplesheet supplies BAM files for it; only prealigned samples keep
// `percent_mapped` in their meta.
//
def readSamplesheet(sample_rows, schema, skip_alignment) {
    validateSamplesheetRows(sample_rows, schema)

    def rows = sample_rows.sort(false) { row -> "${row.fastq_1 ?: row.genome_bam ?: row.transcriptome_bam}" as String }.collect { row ->
        def meta = [ id: row.sample as String, strandedness: row.strandedness ] +
            [ seq_platform: row.seq_platform, seq_center: row.seq_center, percent_mapped: row.percent_mapped ].findAll { _key, value -> value != null }
        def m = meta + [ single_end: !row.fastq_2 ]
        return [ id: m.id, meta: m, run: row.fastq_2 ? [ row.fastq_1, row.fastq_2 ] : [ row.fastq_1 ], genome_bam: row.genome_bam, transcriptome_bam: row.transcriptome_bam ]
    }

    return rows.groupBy { r -> r.id }.collect { id, group ->
        def (meta, runs, genome_bam, transcriptome_bam) = checkSamplesAfterGrouping(
            [ id, group*.meta, group*.run, group*.genome_bam, group*.transcriptome_bam ]
        )
        def prealigned = skip_alignment && (genome_bam || transcriptome_bam) ? true : false
        return [
            id:                id,
            meta:              prealigned ? meta : meta.findAll { key, _value -> key != 'percent_mapped' },
            reads:             prealigned ? [] : runs.flatten(),
            runs:              prealigned ? [] : runs,
            bam:               genome_bam,
            transcriptome_bam: transcriptome_bam,
            percent_mapped:    prealigned ? (meta.percent_mapped ?: null) : null,
            prealigned:        prealigned
        ]
    }
}
