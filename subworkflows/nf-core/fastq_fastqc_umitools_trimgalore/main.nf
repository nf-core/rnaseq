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
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
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
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC.config.ext.args   = params.fastqc_args
        FASTQC (reads)
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    umi_reads = reads
    umi_log   = Channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT.config.ext.args   = [
            params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : '',
            params.umitools_bc_pattern     ? "--bc-pattern='${params.umitools_bc_pattern}'" : '',
            params.umitools_bc_pattern2    ? "--bc-pattern2='${params.umitools_bc_pattern2}'" : '',
            params.umitools_umi_separator  ? "--umi-separator='${params.umitools_umi_separator}'" : ''
        ].join(' ').trim()
        UMITOOLS_EXTRACT.config.publishDir = [
            [
                path: "${params.outdir}/umitools",
                mode: params.publish_dir_mode,
                pattern: "*.log"
            ],
            [
                path: "${params.outdir}/umitools",
                mode: params.publish_dir_mode,
                pattern: "*.fastq.gz",
                enabled: params.save_umi_intermeds
            ]
        ]
        UMITOOLS_EXTRACT (reads)
        umi_reads   = UMITOOLS_EXTRACT.out.reads
        umi_log     = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1,2]) {
            UMITOOLS_EXTRACT
                .out
                .reads
                .map {
                    meta, reads ->
                        meta.single_end ? [ meta, reads ] : [ meta + ['single_end': true], reads[umi_discard_read % 2] ]
                }
                .set { umi_reads }
        }
    }

    trim_reads      = umi_reads
    trim_unpaired   = Channel.empty()
    trim_html       = Channel.empty()
    trim_zip        = Channel.empty()
    trim_log        = Channel.empty()
    trim_read_count = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE.config.ext.args   = {
            [
                "--fastqc_args '-t ${task.cpus}'",
                params.extra_trimgalore_args ? params.extra_trimgalore_args.split("\\s(?=--)") : ''
            ].flatten().unique(false).join(' ').trim()
        }
        TRIMGALORE.config.publishDir = [
            [
                path: "${params.outdir}/${params.trimmer}/fastqc",
                mode: params.publish_dir_mode,
                pattern: "*.{html,zip}"
            ],
            [
                path: "${params.outdir}/${params.trimmer}",
                mode: params.publish_dir_mode,
                pattern: "*.fq.gz",
                enabled: params.save_trimmed
            ],
            [
                path: "${params.outdir}/${params.trimmer}",
                mode: params.publish_dir_mode,
                pattern: "*.txt"
            ]
        ]
        TRIMGALORE (umi_reads)
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_html     = TRIMGALORE.out.html
        trim_zip      = TRIMGALORE.out.zip
        trim_log      = TRIMGALORE.out.log
        ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        TRIMGALORE
            .out
            .reads
            .join(trim_log, remainder: true)
            .map {
                meta, reads, trim_log ->
                    if (trim_log) {
                        num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? trim_log : trim_log[-1])
                        [ meta, reads, num_reads ]
                    } else {
                        [ meta, reads, min_trimmed_reads.toFloat() + 1 ]
                    }
            }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { meta, reads, num_reads -> num_reads >= min_trimmed_reads.toFloat() }
            .map { meta, reads, num_reads -> [ meta, reads ] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, reads, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
