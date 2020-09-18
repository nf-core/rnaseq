/*
 * Read QC, UMI extraction and trimming
 */

include { FASTQC           } from '../software/fastqc/main'
include { UMITOOLS_EXTRACT } from '../software/umitools/extract/main'
include { TRIMGALORE       } from '../software/trimgalore/main'

workflow FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    skip_fastqc        // boolean: true/false
    umi_extract        // boolean: true/false
    skip_trimming      // boolean: true/false

    fastqc_options     // map: options for FastQC module
    umitools_options   // map: options for UMI-tools extract module
    trimgalore_options // map: options for TrimGalore! module

    main:
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC ( reads, fastqc_options ).html.set { fastqc_html }
        fastqc_zip = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    umi_reads = reads
    umi_log = Channel.empty()
    umitools_version = Channel.empty()
    if (umi_extract) {
        UMITOOLS_EXTRACT ( reads, umitools_options ).reads.set { umi_reads }
        umi_log = UMITOOLS_EXTRACT.out.log
        umitools_version = UMITOOLS_EXTRACT.out.version
    }

    trim_reads = umi_reads
    trim_html = Channel.empty()
    trim_zip = Channel.empty()
    trim_log = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE ( umi_reads, trimgalore_options ).reads.set { trim_reads }
        trim_html = TRIMGALORE.out.html
        trim_zip = TRIMGALORE.out.zip
        trim_log = TRIMGALORE.out.log
        trimgalore_version = TRIMGALORE.out.version
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]
    fastqc_version     //    path: *.version.txt

    umi_log            // channel: [ val(meta), [ log ] ]
    umitools_version   //    path: *.version.txt

    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trimgalore_version //    path: *.version.txt
}
