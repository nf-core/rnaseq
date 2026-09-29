//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'
include { FastqFastqcUmitoolsFastp } from './types'

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
        fastqc_raw_html = FASTQC_RAW.out.map { r -> [r.meta, r.html] }
        fastqc_raw_zip = FASTQC_RAW.out.map { r -> [r.meta, r.zip] }

        ch_results = ch_results
            .join(FASTQC_RAW.out.map { r -> [r.id, r] }, by: [0])
            .map { id, fields, r ->
                [id, fields + [fastqc: record(raw_html: r.html, raw_zip: r.zip, trim_html: null, trim_zip: null)]]
            }
    }

    trimmer_reads = reads_only
    umi_reads = channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(
            reads_only
        )
        // Typed processes read a lone Path as its name components, so wrap single-end reads in a list
        trimmer_reads = UMITOOLS_EXTRACT.out.reads.map { meta, reads_ -> [meta, [reads_].flatten()] }
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
                    meta.single_end ? [meta, [_reads].flatten()] : [meta + [single_end: true], [_reads[umi_discard_read % 2]]]
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
                reads.map { meta, _original_reads, adapter_fasta -> [meta.id, adapter_fasta ?: null] }
            )
            .map { _sample_id, meta, umi_reads_files, adapter_fasta -> [meta, umi_reads_files, adapter_fasta] }

        FASTP(
            umi_reads_with_adapters,
            false,
            save_trimmed_fail,
            save_merged
        )
        trim_json = FASTP.out.map { r -> [r.meta, r.json] }
        trim_html = FASTP.out.map { r -> [r.meta, r.html] }
        trim_log = FASTP.out.map { r -> [r.meta, r.log] }
        trim_reads_fail = FASTP.out.filter { r -> r.reads_fail }.map { r -> [r.meta, r.reads_fail] }
        trim_reads_merged = FASTP.out.filter { r -> r.reads_merged }.map { r -> [r.meta, r.reads_merged] }

        // FASTP reads are optional, so a sample can lack a read count
        FASTP.out
            .map { r ->
                def num_reads = r.reads ? getFastpReadsAfterFiltering(r.json, min_trimmed_reads.toLong()) : null
                [r, num_reads, getFastpAdapterSequence(r.json)]
            }
            .set { ch_trimmed }

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_trimmed
            .filter { r, num_reads, _seq -> r.reads && num_reads >= min_trimmed_reads.toLong() }
            .map { r, _num_reads, _seq -> [r.meta, r.reads] }
            .set { trim_reads }

        ch_trimmed
            .filter { r, _num_reads, _seq -> r.reads }
            .map { r, num_reads, _seq -> [r.meta, num_reads] }
            .set { trim_read_count }

        ch_trimmed
            .map { r, _num_reads, seq -> [r.meta, seq] }
            .set { adapter_seq }

        ch_results = ch_results
            .join(ch_trimmed.map { r, num_reads, seq -> [r.id, r, num_reads, seq] }, by: [0])
            .map { id, fields, r, num_reads, seq ->
                def trim = record(
                    html:         r.html,
                    json:         r.json,
                    log:          r.log,
                    reads_fail:   r.reads_fail ?: null,
                    reads_merged: r.reads_merged
                )
                [id, fields + [trim: trim, adapter_seq: seq, num_trimmed_reads: num_reads != null ? num_reads as Long : null]]
            }

        if (!skip_fastqc) {
            FASTQC_TRIM(
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.map { r -> [r.meta, r.html] }
            fastqc_trim_zip = FASTQC_TRIM.out.map { r -> [r.meta, r.zip] }

            // Samples below min_trimmed_reads are not passed to FASTQC_TRIM
            ch_results = ch_results
                .join(FASTQC_TRIM.out.map { r -> [r.id, r] }, by: [0], remainder: true)
                .map { id, fields, r ->
                    def fastqc = record(
                        raw_html:  fields.fastqc.raw_html,
                        raw_zip:   fields.fastqc.raw_zip,
                        trim_html: r ? r.html : null,
                        trim_zip:  r ? r.zip : null
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
        )
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
