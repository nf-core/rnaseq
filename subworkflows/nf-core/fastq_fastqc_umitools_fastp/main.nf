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
    catch (Exception ex) {
        return ""
    }
}

workflow FASTQ_FASTQC_UMITOOLS_FASTP {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    skip_trimming     // boolean: true/false
    adapter_fasta     // file: adapter.fasta
    save_trimmed_fail // boolean: true/false
    save_merged       // boolean: true/false
    min_trimmed_reads // integer: > 0

    main:
    ch_versions = Channel.empty()
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip = Channel.empty()
    umi_log = Channel.empty()
    trim_json = Channel.empty()
    trim_html = Channel.empty()
    trim_log = Channel.empty()
    trim_reads_fail = Channel.empty()
    trim_reads_merged = Channel.empty()
    fastqc_trim_html = Channel.empty()
    fastqc_trim_zip = Channel.empty()
    trim_read_count = Channel.empty()
    adapter_seq = Channel.empty()

    if (!skip_fastqc) {
        FASTQC_RAW(
            reads
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip = FASTQC_RAW.out.zip
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    umi_reads = reads
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(
            reads
        )
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, _reads ->
                    meta.single_end ? [meta, _reads] : [meta + [single_end: true], _reads[umi_discard_read % 2]]
                }
                .set { umi_reads }
        }
    }

    trim_reads = umi_reads
    if (!skip_trimming) {
        FASTP(
            umi_reads,
            adapter_fasta,
            false,
            save_trimmed_fail,
            save_merged,
        )
        trim_json = FASTP.out.json
        trim_html = FASTP.out.html
        trim_log = FASTP.out.log
        trim_reads_fail = FASTP.out.reads_fail
        trim_reads_merged = FASTP.out.reads_merged
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

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
            ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        }
    }

    emit:
    reads             = trim_reads // channel: [ val(meta), [ reads ] ]
    fastqc_raw_html   // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip    // channel: [ val(meta), [ zip ] ]
    umi_log           // channel: [ val(meta), [ log ] ]
    adapter_seq       // channel: [ val(meta), [ adapter_seq] ]
    trim_json         // channel: [ val(meta), [ json ] ]
    trim_html         // channel: [ val(meta), [ html ] ]
    trim_log          // channel: [ val(meta), [ log ] ]
    trim_reads_fail   // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged // channel: [ val(meta), [ fastq.gz ] ]
    trim_read_count   // channel: [ val(meta), val(count) ]
    fastqc_trim_html  // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip   // channel: [ val(meta), [ zip ] ]
    versions          = ch_versions // channel: [ versions.yml ]
}
