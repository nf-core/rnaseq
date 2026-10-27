include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'


workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams    // channel: [ meta, srbam, lrbam ]
    ch_mode    // channel: [ val(mode) ]
    ch_chrgtf  // channel: [ meta, gtf ]

    main:

    STRINGTIE_STRINGTIE(
        ch_bams,
        ch_mode,
        ch_chrgtf.map { _meta, gtf -> [gtf] }
    )

    STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { _meta, gtf -> gtf }
        .toSortedList { a, b -> a.name <=> b.name }
        .filter { gtfs -> gtfs.size() > 0 }
        .map { gtfs -> [ [id: 'stringtie_merge'], gtfs ] }
        .set { collected_gtfs }

    STRINGTIE_MERGE(
        collected_gtfs,
        ch_chrgtf
    )

    emit:
    stringtie_gtf = STRINGTIE_MERGE.out.merged_gtf // channel: [ meta, gtf ]
}
