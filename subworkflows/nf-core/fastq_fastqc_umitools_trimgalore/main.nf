nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'
include { TRIMGALORE       } from '../../../modules/nf-core/trimgalore/main'
include { ReadsInput; FastqcResult; FastqFastqcUmitoolsTrimgalore } from '../../../modules/nf-core/types'

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s*([\d\.]+)/
        if (total_reads_matcher) {
            total_reads = total_reads_matcher[0][1].toFloat()
        }
        if (filtered_reads_matcher) {
            filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
    }
    return total_reads - filtered_reads
}

workflow FASTQ_FASTQC_UMITOOLS_TRIMGALORE {
    take:
    ch_reads: Channel<ReadsInput>
    skip_fastqc: Boolean // true/false
    with_umi: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    skip_trimming: Boolean // true/false
    umi_discard_read: Integer // 0, 1 or 2
    min_trimmed_reads: Integer // > 0

    main:
    // Each stage that runs joins its outputs onto this per-sample record, overwriting
    // the null placeholders of the fields it owns.
    def ch_results: Channel<FastqFastqcUmitoolsTrimgalore> = ch_reads.map { r ->
        record(
            id:                r.id,
            meta:              r.meta,
            reads:             r.reads,
            fastqc_raw_html:   null,
            fastqc_raw_zip:    null,
            umi:               null,
            trim:              null,
            num_trimmed_reads: null
        )
    }

    if (!skip_fastqc) {
        def ch_fastqc: Channel<FastqcResult> = FASTQC(ch_reads)
        ch_results = ch_results.join(
            ch_fastqc.map { r -> record(id: r.id, fastqc_raw_html: r.html, fastqc_raw_zip: r.zip) },
            by: 'id'
        )
    }

    ch_trimmer_reads = ch_reads.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, umi: null) }
    if (with_umi && !skip_umi_extract) {
        ch_trimmer_reads = UMITOOLS_EXTRACT(ch_reads).map { r ->
            // Discard R1 / R2 if required
            def discard = umi_discard_read in [1, 2] && !r.meta.single_end
            def meta = r.meta
            if (discard) {
                meta = r.meta + [single_end: true]
            }
            record(
                id:    r.id,
                meta:  meta,
                reads: discard ? [r.reads[umi_discard_read % 2]] : r.reads,
                umi:   record(log: r.log, reads: r.reads)
            )
        }
    }
    ch_results = ch_results.join(ch_trimmer_reads, by: 'id')

    if (!skip_trimming) {
        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_trim = TRIMGALORE(ch_trimmer_reads).map { r ->
            def num_reads = r.log.isEmpty()
                ? (min_trimmed_reads as Float) + 1
                : getTrimGaloreReadsAfterFiltering(r.log[-1]) as Float
            record(
                id:                r.id,
                meta:              r.meta,
                reads:             num_reads >= (min_trimmed_reads as Float) ? r.reads : null,
                // TrimGalore reports and unpaired reads are all optional outputs
                trim: record(
                    html:     r.html.isEmpty() ? null : r.html,
                    zip:      r.zip.isEmpty() ? null : r.zip,
                    log:      r.log.isEmpty() ? null : r.log,
                    json:     r.json.isEmpty() ? null : r.json,
                    unpaired: r.unpaired.isEmpty() ? null : r.unpaired
                ),
                num_trimmed_reads: num_reads
            )
        }
        ch_results = ch_results.join(ch_trim, by: 'id')
    }

    emit:
    ch_results
}
