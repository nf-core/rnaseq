//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../modules/nf-core/modules/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/modules/umitools/extract/main'
include { TRIMGALORE       } from '../../modules/nf-core/modules/trimgalore/main'

workflow FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    with_umi      // boolean: true/false
    skip_trimming // boolean: true/false

    main:

    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC ( reads ).html.set { fastqc_html }
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    umi_reads = reads
    umi_log   = Channel.empty()
    if (with_umi) {
        UMITOOLS_EXTRACT ( reads ).reads.set { umi_reads }
        umi_log     = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
    }

    trim_reads = umi_reads
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE ( umi_reads ).reads.set { trim_reads }
        trim_html   = TRIMGALORE.out.html
        trim_zip    = TRIMGALORE.out.zip
        trim_log    = TRIMGALORE.out.log
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
