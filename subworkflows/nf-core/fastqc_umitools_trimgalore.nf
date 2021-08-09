//
// Read QC, UMI extraction and trimming
//

params.fastqc_options     = [:]
params.umitools_options   = [:]
params.trimgalore_options = [:]

include { FASTQC           } from '../../modules/nf-core/modules/fastqc/main'           addParams( options: params.fastqc_options     )
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/modules/umitools/extract/main' addParams( options: params.umitools_options   )
include { TRIMGALORE       } from '../../modules/nf-core/modules/trimgalore/main'       addParams( options: params.trimgalore_options )

workflow FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    with_umi      // boolean: true/false
    skip_trimming // boolean: true/false

    main:
    fastqc_html    = Channel.empty()
    fastqc_zip     = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC ( reads ).html.set { fastqc_html }
        fastqc_zip     = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    umi_reads        = reads
    umi_log          = Channel.empty()
    umitools_version = Channel.empty()
    if (with_umi) {
        UMITOOLS_EXTRACT ( reads ).reads.set { umi_reads }
        umi_log          = UMITOOLS_EXTRACT.out.log
        umitools_version = UMITOOLS_EXTRACT.out.version
    }

    trim_reads = umi_reads
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE ( umi_reads ).reads.set { trim_reads }
        trim_html  = TRIMGALORE.out.html
        trim_zip   = TRIMGALORE.out.zip
        trim_log   = TRIMGALORE.out.log
        trimgalore_version = TRIMGALORE.out.version
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]
    fastqc_version     //    path: versions.yml

    umi_log            // channel: [ val(meta), [ log ] ]
    umitools_version   //    path: versions.yml

    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trimgalore_version //    path: versions.yml
}
