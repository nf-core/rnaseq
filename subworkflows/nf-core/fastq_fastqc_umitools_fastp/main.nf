nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'
include { FastpReads; FastqFastqcUmitoolsFastp } from '../../../modules/nf-core/types'
include { FastqcResult } from '../../../modules/nf-core/fastqc/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//

def getFastpReadsAfterFiltering(json_file, min_num_reads) {

    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = (new groovy.json.JsonSlurper().parseText(json_file.text) as Map)['summary'] as Map
    return (json['after_filtering'] as Map)['total_reads'] as Long
}

def getFastpAdapterSequence(json_file) {
    // Handle stub runs
    if (workflow.stubRun) {
        return ""
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text) as Map
    try {
        return (json['adapter_cutting'] as Map)['read1_adapter_sequence'] as String
    }
    catch (Exception _ex) {
        return ""
    }
}

workflow FASTQ_FASTQC_UMITOOLS_FASTP {
    take:
    ch_reads: Channel<FastpReads>
    skip_fastqc: Boolean // true/false
    with_umi: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    umi_discard_read: Integer // 0, 1 or 2
    skip_trimming: Boolean // true/false
    save_trimmed_fail: Boolean // true/false
    save_merged: Boolean // true/false
    min_trimmed_reads: Integer // > 0

    main:
    // Each stage that runs joins its outputs onto this per-sample record, overwriting
    // the null placeholders of the fields it owns.
    def ch_results: Channel<FastqFastqcUmitoolsFastp> = ch_reads.map { r ->
        record(
            id:                r.id,
            meta:              r.meta,
            reads:             r.reads,
            fastqc_raw_html:   null,
            fastqc_raw_zip:    null,
            fastqc_trim_html:  null,
            fastqc_trim_zip:   null,
            umi:               null,
            trim:              null,
            adapter_seq:       null,
            num_trimmed_reads: null
        )
    }

    if (!skip_fastqc) {
        def ch_fastqc_raw: Channel<FastqcResult> = FASTQC_RAW(ch_reads)
        ch_results = ch_results.join(
            ch_fastqc_raw.map { r -> record(id: r.id, fastqc_raw_html: r.html, fastqc_raw_zip: r.zip) },
            by: 'id'
        )
    }

    ch_trimmer_reads = ch_reads.map { r ->
        record(id: r.id, meta: r.meta, reads: r.reads, adapter_fasta: r.adapter_fasta, umi: null)
    }
    if (with_umi && !skip_umi_extract) {
        // The adapter fasta of the original input is re-attached by sample id, since UMI extraction does not carry it
        ch_trimmer_reads = UMITOOLS_EXTRACT(ch_reads)
            .join(ch_reads.map { r -> record(id: r.id, adapter_fasta: r.adapter_fasta) }, by: 'id')
            .map { r ->
                // Discard R1 / R2 if required
                def discard = umi_discard_read in [1, 2] && !r.meta.single_end
                def meta = r.meta
                if (discard) {
                    meta = r.meta + [single_end: true]
                }
                record(
                    id:            r.id,
                    meta:          meta,
                    reads:         discard ? [r.reads[umi_discard_read % 2]] : r.reads,
                    adapter_fasta: r.adapter_fasta,
                    umi:           record(log: r.log, reads: r.reads)
                )
            }
    }
    ch_results = ch_results.join(
        ch_trimmer_reads.map { r -> record(id: r.id, meta: r.meta, reads: r.reads, umi: r.umi) },
        by: 'id'
    )

    if (!skip_trimming) {
        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_trim = FASTP(ch_trimmer_reads, false, save_trimmed_fail, save_merged).map { r ->
            // FASTP reads are optional, so a sample can lack a read count
            def num_reads = r.reads.isEmpty() ? null : getFastpReadsAfterFiltering(r.json, min_trimmed_reads as Long)
            record(
                id:                r.id,
                meta:              r.meta,
                reads:             num_reads != null && num_reads >= (min_trimmed_reads as Long) ? r.reads : null,
                trim:              record(
                    html:         r.html,
                    json:         r.json,
                    log:          r.log,
                    reads_fail:   r.reads_fail.isEmpty() ? null : r.reads_fail,
                    reads_merged: r.reads_merged
                ),
                adapter_seq:       getFastpAdapterSequence(r.json),
                num_trimmed_reads: num_reads
            )
        }
        ch_results = ch_results.join(ch_trim, by: 'id')

        if (!skip_fastqc) {
            def ch_fastqc_trim: Channel<FastqcResult> = FASTQC_TRIM(ch_results.filter { r -> r.reads != null })

            // Samples below min_trimmed_reads are not passed to FASTQC_TRIM
            ch_results = ch_results.join(
                ch_fastqc_trim.map { r -> record(id: r.id, fastqc_trim_html: r.html, fastqc_trim_zip: r.zip) },
                by: 'id',
                remainder: true
            )
        }
    }

    emit:
    ch_results
}
