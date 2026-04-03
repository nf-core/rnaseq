include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE    } from '../../../modules/nf-core/stringtie/merge/main'


workflow BAM_STRINGTIE_MERGE {

    take:
    bam_sorted // channel: [ meta, bam ]
    ch_chrgtf  // channel: [ meta, gtf ]

    main:
    ch_versions       = channel.empty()
    ch_stringtie_gtfs = channel.empty()

    STRINGTIE_STRINGTIE(
        bam_sorted,
        ch_chrgtf.map { _meta, gtf -> [ gtf ] }
    )

    STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { _meta, gtf -> [ gtf ] }
        .set { stringtie_gtfs }

    STRINGTIE_MERGE(
        stringtie_gtfs,
        ch_chrgtf.map { _meta, gtf -> [ gtf ] }
    )
    ch_stringtie_gtfs = STRINGTIE_MERGE.out.gtf

    emit:
    stringtie_gtf = ch_stringtie_gtfs // channel: [ meta, gtf ]
    versions      = ch_versions                     // channel: [ path(versions.yml) ]
}
