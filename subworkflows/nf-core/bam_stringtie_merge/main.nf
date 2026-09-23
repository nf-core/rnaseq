include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'

record StringtieAssembly {
    id:             String
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
}

record StringtieMerged {
    id:         String
    merged_gtf: Path
}

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

    ch_results = STRINGTIE_STRINGTIE.out.transcript_gtf
        .join(STRINGTIE_STRINGTIE.out.abundance)
        .join(STRINGTIE_STRINGTIE.out.coverage_gtf, remainder: true)
        .join(STRINGTIE_STRINGTIE.out.ballgown, remainder: true)
        .map { meta, transcript_gtf, abundance, coverage_gtf, ballgown ->
            record(
                id:             meta.id,
                transcript_gtf: transcript_gtf,
                abundance:      abundance,
                coverage_gtf:   coverage_gtf,
                ballgown:       ballgown ? [ballgown].flatten() : null
            ) as StringtieAssembly
        }

    ch_merged_results = STRINGTIE_MERGE.out.merged_gtf
        .map { meta, gtf -> record(id: meta.id, merged_gtf: gtf) as StringtieMerged }

    emit:
    stringtie_gtf  = STRINGTIE_MERGE.out.merged_gtf // channel: [ meta, gtf ]
    transcript_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf // channel: [ meta, gtf ]
    abundance      = STRINGTIE_STRINGTIE.out.abundance // channel: [ meta, txt ]
    coverage_gtf   = STRINGTIE_STRINGTIE.out.coverage_gtf // channel: [ meta, gtf ]
    ballgown       = STRINGTIE_STRINGTIE.out.ballgown // channel: [ meta, [ ctab ] ]
    results        = ch_results // channel: StringtieAssembly
    merged_results = ch_merged_results // channel: StringtieMerged
}
