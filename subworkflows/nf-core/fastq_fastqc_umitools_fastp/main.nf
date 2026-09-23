//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//

def getFastpReadsAfterFiltering(json_file, min_num_reads) {

    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary') as Map
    return json['after_filtering']['total_reads'].toLong()
}

def getFastpAdapterSequence(json_file) {
    // Handle stub runs
    if (workflow.stubRun) {
        return ""
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text) as Map
    try {
        return json['adapter_cutting']['read1_adapter_sequence']
    }
    catch (Exception _ex) {
        return ""
    }
}

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

workflow FASTQ_FASTQC_UMITOOLS_FASTP {
    take:
    reads             // channel: [ val(meta), [ reads ], adapter_fasta ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    skip_trimming     // boolean: true/false
    save_trimmed_fail // boolean: true/false
    save_merged       // boolean: true/false
    min_trimmed_reads // integer: > 0

    main:
    fastqc_raw_html = channel.empty()
    fastqc_raw_zip = channel.empty()
    umi_log = channel.empty()
    trim_json = channel.empty()
    trim_html = channel.empty()
    trim_log = channel.empty()
    trim_reads_fail = channel.empty()
    trim_reads_merged = channel.empty()
    fastqc_trim_html = channel.empty()
    fastqc_trim_zip = channel.empty()
    trim_read_count = channel.empty()
    adapter_seq = channel.empty()

    // Split input channel for reads-only operations
    reads_only = reads.map { meta, reads_files, _adapter_fasta -> [ meta, reads_files ] }

    // Each step that runs joins its outputs onto this per-sample skeleton;
    // groups for skipped steps are left unset and become null in the final
    // record, so no join is ever made against an empty channel.
    ch_results = reads_only.map { meta, _reads -> [meta.id, [:]] }

    if (!skip_fastqc) {
        FASTQC_RAW(
            reads_only
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip = FASTQC_RAW.out.zip

        ch_results = ch_results
            .join(FASTQC_RAW.out.html.map { meta, html -> [meta.id, html] }, by: [0])
            .join(FASTQC_RAW.out.zip.map { meta, zip -> [meta.id, zip] }, by: [0])
            .map { id, fields, html, zip ->
                [id, fields + [fastqc: record(raw_html: [html].flatten(), raw_zip: [zip].flatten(), trim_html: null, trim_zip: null)]]
            }
    }

    trimmer_reads = reads_only
    umi_reads = channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(
            reads_only
        )
        trimmer_reads = UMITOOLS_EXTRACT.out.reads
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log

        ch_results = ch_results
            .join(UMITOOLS_EXTRACT.out.log.map { meta, log -> [meta.id, log] }, by: [0])
            .join(UMITOOLS_EXTRACT.out.reads.map { meta, reads_ -> [meta.id, reads_] }, by: [0])
            .map { id, fields, log, reads_ ->
                [id, fields + [umi: record(log: log, reads: [reads_].flatten())]]
            }

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, _reads ->
                    meta.single_end ? [meta, _reads] : [meta + [single_end: true], _reads[umi_discard_read % 2]]
                }
                .set { trimmer_reads }
        }
    }

    trim_reads = trimmer_reads
    if (!skip_trimming) {
        // Rejoin trimmer_reads with adapter info from original input
        // Use ID-based join to handle metadata modifications from UMI processing
        umi_reads_with_adapters = trimmer_reads
            .map { meta, reads_files -> [meta.id, meta, reads_files] }
            .join(
                reads.map { meta, _original_reads, adapter_fasta -> [meta.id, adapter_fasta ?: []] }
            )
            .map { _sample_id, meta, umi_reads_files, adapter_fasta -> [meta, umi_reads_files, adapter_fasta] }

        FASTP(
            umi_reads_with_adapters,
            false,
            save_trimmed_fail,
            save_merged
        )
        trim_json = FASTP.out.json
        trim_html = FASTP.out.html
        trim_log = FASTP.out.log
        trim_reads_fail = FASTP.out.reads_fail
        trim_reads_merged = FASTP.out.reads_merged

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        FASTP.out.reads.join(trim_json).map { meta, _reads, json -> [meta, _reads, getFastpReadsAfterFiltering(json, min_trimmed_reads.toLong())] }.set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { _meta, _reads, num_reads -> num_reads >= min_trimmed_reads.toLong() }
            .map { meta, _reads, _num_reads -> [meta, _reads] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [meta, num_reads] }
            .set { trim_read_count }

        trim_json
            .map { meta, json -> [meta, getFastpAdapterSequence(json)] }
            .set { adapter_seq }

        // FASTP.out.reads is optional, so a sample can lack a read count
        ch_results = ch_results
            .join(trim_json.map { meta, f -> [meta.id, f] }, by: [0])
            .join(trim_html.map { meta, f -> [meta.id, f] }, by: [0])
            .join(trim_log.map { meta, f -> [meta.id, f] }, by: [0])
            .join(adapter_seq.map { meta, seq -> [meta.id, seq] }, by: [0])
            .join(trim_reads_fail.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(trim_reads_merged.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(trim_read_count.map { meta, num_reads -> [meta.id, num_reads] }, by: [0], remainder: true)
            .map { id, fields, json, html, log, seq, reads_fail, reads_merged, num_reads ->
                def trim = record(
                    html:         html,
                    json:         json,
                    log:          log,
                    reads_fail:   reads_fail ? [reads_fail].flatten() : null,
                    reads_merged: reads_merged
                )
                [id, fields + [trim: trim, adapter_seq: seq, num_trimmed_reads: num_reads != null ? num_reads as Long : null]]
            }

        if (!skip_fastqc) {
            FASTQC_TRIM(
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip = FASTQC_TRIM.out.zip

            // Samples below min_trimmed_reads are not passed to FASTQC_TRIM
            ch_results = ch_results
                .join(FASTQC_TRIM.out.html.map { meta, html -> [meta.id, html] }, by: [0], remainder: true)
                .join(FASTQC_TRIM.out.zip.map { meta, zip -> [meta.id, zip] }, by: [0], remainder: true)
                .map { id, fields, html, zip ->
                    def fastqc = record(
                        raw_html:  fields.fastqc.raw_html,
                        raw_zip:   fields.fastqc.raw_zip,
                        trim_html: html ? [html].flatten() : null,
                        trim_zip:  zip ? [zip].flatten() : null
                    )
                    [id, fields + [fastqc: fastqc]]
                }
        }
    }

    ch_results = ch_results.map { id, fields ->
        record(
            id:                id,
            fastqc:            fields.fastqc,
            umi:               fields.umi,
            trim:              fields.trim,
            adapter_seq:       fields.adapter_seq,
            num_trimmed_reads: fields.num_trimmed_reads
        ) as FastqFastqcUmitoolsFastp
    }

    emit:
    reads             = trim_reads // channel: [ val(meta), [ reads ] ]
    fastqc_raw_html   // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip    // channel: [ val(meta), [ zip ] ]
    umi_log           // channel: [ val(meta), [ log ] ]
    umi_reads         // channel: [ val(meta), [ reads ] ]
    adapter_seq       // channel: [ val(meta), [ adapter_seq] ]
    trim_json         // channel: [ val(meta), [ json ] ]
    trim_html         // channel: [ val(meta), [ html ] ]
    trim_log          // channel: [ val(meta), [ log ] ]
    trim_reads_fail   // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged // channel: [ val(meta), [ fastq.gz ] ]
    trim_read_count   // channel: [ val(meta), val(count) ]
    fastqc_trim_html  // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip   // channel: [ val(meta), [ zip ] ]
    results           = ch_results // channel: FastqFastqcUmitoolsFastp
}
