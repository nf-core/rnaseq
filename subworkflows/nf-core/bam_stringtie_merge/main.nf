include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'
include { StringtieAssembly; StringtieMerged } from './types'

workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams    // channel: [ meta, srbam, lrbam ]
    ch_mode    // channel: [ val(mode) ], a list of mode names
    ch_chrgtf  // channel: [ meta, gtf ]

    main:

    STRINGTIE_STRINGTIE(
        ch_bams,
        ch_mode,
        ch_chrgtf.map { _meta, gtf -> gtf }
    )

    STRINGTIE_STRINGTIE.out
        .map { r -> r.transcript_gtf }
        .toSortedList { a, b -> a.name <=> b.name }
        .filter { gtfs -> gtfs.size() > 0 }
        .map { gtfs -> [ [id: 'stringtie_merge'], gtfs ] }
        .set { collected_gtfs }

    STRINGTIE_MERGE(
        collected_gtfs,
        ch_chrgtf
    )

    emit:
    results        = STRINGTIE_STRINGTIE.out // channel: StringtieAssembly
    merged_results = STRINGTIE_MERGE.out     // channel: StringtieMerged
}
