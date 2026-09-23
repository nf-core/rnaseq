//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'
include { TRIMGALORE       } from '../../../modules/nf-core/trimgalore/main'

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

workflow FASTQ_FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    skip_trimming     // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    min_trimmed_reads // integer: > 0

    main:
    // Each step that runs joins its outputs onto this per-sample skeleton;
    // groups for skipped steps are left unset and become null in the final
    // record, so no join is ever made against an empty channel.
    ch_results = reads.map { meta, _reads -> [meta.id, [:]] }

    ch_fastqc_html = channel.empty()
    ch_fastqc_zip = channel.empty()
    if (!skip_fastqc) {
        FASTQC(reads)
        ch_fastqc_html = FASTQC.out.html
        ch_fastqc_zip = FASTQC.out.zip

        ch_results = ch_results
            .join(FASTQC.out.html.map { meta, html -> [meta.id, html] }, by: [0])
            .join(FASTQC.out.zip.map { meta, zip -> [meta.id, zip] }, by: [0])
            .map { id, fields, html, zip ->
                [id, fields + [fastqc: record(raw_html: [html].flatten(), raw_zip: [zip].flatten())]]
            }
    }

    ch_trimmer_reads = reads
    ch_umi_log = channel.empty()
    ch_umi_reads = channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(reads)
        ch_trimmer_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log

        ch_results = ch_results
            .join(UMITOOLS_EXTRACT.out.log.map { meta, log -> [meta.id, log] }, by: [0])
            .join(UMITOOLS_EXTRACT.out.reads.map { meta, reads_ -> [meta.id, reads_] }, by: [0])
            .map { id, fields, log, reads_ ->
                [id, fields + [umi: record(log: log, reads: [reads_].flatten())]]
            }

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, reads_ ->
                    meta.single_end ? [meta, reads_] : [meta + ['single_end': true], reads_[umi_discard_read % 2]]
                }
                .set { ch_trimmer_reads }
        }
    }

    ch_trim_reads = ch_trimmer_reads
    ch_trim_unpaired = channel.empty()
    ch_trim_html = channel.empty()
    ch_trim_zip = channel.empty()
    ch_trim_log = channel.empty()
    ch_trim_json = channel.empty()
    ch_trim_read_count = channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(ch_trimmer_reads)
        ch_trim_unpaired = TRIMGALORE.out.unpaired
        ch_trim_html = TRIMGALORE.out.html
        ch_trim_zip = TRIMGALORE.out.zip
        ch_trim_log = TRIMGALORE.out.log
        ch_trim_json = TRIMGALORE.out.json

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        TRIMGALORE.out.reads
            .join(ch_trim_log, remainder: true)
            .map { meta, reads_, trim_log ->
                if (trim_log) {
                    def num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? trim_log : trim_log[-1])
                    [meta, reads_, num_reads]
                }
                else {
                    [meta, reads_, min_trimmed_reads.toFloat() + 1]
                }
            }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { _meta, _reads, num_reads -> num_reads >= min_trimmed_reads.toFloat() }
            .map { meta, reads_, _num_reads -> [meta, reads_] }
            .set { ch_trim_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [meta, num_reads] }
            .set { ch_trim_read_count }

        // TrimGalore reports and unpaired reads are all optional module outputs
        ch_results = ch_results
            .join(ch_trim_read_count.map { meta, num_reads -> [meta.id, num_reads] }, by: [0])
            .join(ch_trim_html.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(ch_trim_zip.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(ch_trim_log.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(ch_trim_json.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(ch_trim_unpaired.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .map { id, fields, num_reads, html, zip, log, json, unpaired ->
                def trim = record(
                    html:     html ? [html].flatten() : null,
                    zip:      zip ? [zip].flatten() : null,
                    log:      log ? [log].flatten() : null,
                    json:     json ? [json].flatten() : null,
                    unpaired: unpaired ? [unpaired].flatten() : null
                )
                [id, fields + [trim: trim, num_trimmed_reads: num_reads as Float]]
            }
    }

    ch_results = ch_results.map { id, fields ->
        record(
            id:                id,
            fastqc:            fields.fastqc,
            umi:               fields.umi,
            trim:              fields.trim,
            num_trimmed_reads: fields.num_trimmed_reads
        ) as FastqFastqcUmitoolsTrimgalore
    }

    emit:
    reads           = ch_trim_reads       // channel: [ val(meta), [ reads ] ]
    fastqc_html     = ch_fastqc_html      // channel: [ val(meta), [ html ] ]
    fastqc_zip      = ch_fastqc_zip       // channel: [ val(meta), [ zip ] ]
    umi_log         = ch_umi_log          // channel: [ val(meta), [ log ] ]
    umi_reads       = ch_umi_reads        // channel: [ val(meta), [ reads ] ]
    trim_unpaired   = ch_trim_unpaired    // channel: [ val(meta), [ reads ] ]
    trim_html       = ch_trim_html        // channel: [ val(meta), [ html ] ]
    trim_zip        = ch_trim_zip         // channel: [ val(meta), [ zip ] ]
    trim_log        = ch_trim_log         // channel: [ val(meta), [ txt ] ]
    trim_json       = ch_trim_json        // channel: [ val(meta), [ json ] ]
    trim_read_count = ch_trim_read_count  // channel: [ val(meta), val(count) ]
    results         = ch_results          // channel: FastqFastqcUmitoolsTrimgalore
}
