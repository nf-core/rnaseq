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

    if (!skip_fastqc) {
        FASTQC_RAW(
            reads_only
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip = FASTQC_RAW.out.zip
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

        if (!skip_fastqc) {
            FASTQC_TRIM(
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip = FASTQC_TRIM.out.zip
        }
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
}
