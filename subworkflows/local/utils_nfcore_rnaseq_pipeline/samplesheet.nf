//
// Samplesheet reading. Kept untyped because `nextflow lint` rejects calls to plugin functions such
// as `samplesheetToList` in typed scripts (nextflow-io/nextflow#7720).
//

include { samplesheetToList         } from 'plugin/nf-schema'
include { checkSamplesAfterGrouping } from './main'

//
// Validate the samplesheet against its schema and return its rows, one per sequencing run:
// [ meta, fastq_1, fastq_2, genome_bam, transcriptome_bam ].
//
def loadSamplesheet(samplesheet, schema) {
    return samplesheetToList(samplesheet, schema)
}

//
// Group the samplesheet rows into one plain map per sample, merging the rows of a sample that has several
// sequencing runs. Keys: id, meta, reads (all FASTQ files), runs (one list of FASTQ files per run),
// bam, transcriptome_bam, percent_mapped and prealigned. A sample is prealigned when alignment is
// skipped and the samplesheet supplies BAM files for it; only prealigned samples keep
// `percent_mapped` in their meta.
//
def readSamplesheet(samplesheet_rows, skip_alignment) {
    def rows = samplesheet_rows.collect { row ->
        def (meta, fastq_1, fastq_2, genome_bam, transcriptome_bam) = row
        def m = meta + [ id: meta.id as String, single_end: !fastq_2 ]
        return [ id: m.id, meta: m, run: fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ], genome_bam: genome_bam, transcriptome_bam: transcriptome_bam ]
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
