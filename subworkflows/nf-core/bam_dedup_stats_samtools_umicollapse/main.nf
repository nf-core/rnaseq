//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE    } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai/csi) ]

    main:

    ch_versions = Channel.empty()

    //
    // umicollapse in bam mode (thus hardcode mode input channel to 'bam')
    //
    UMICOLLAPSE ( ch_bam_bai, channel.value( 'bam' ))
    ch_versions = ch_versions.mix(UMICOLLAPSE.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( UMICOLLAPSE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_bai_dedup = UMICOLLAPSE.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }

    BAM_STATS_SAMTOOLS ( ch_bam_bai_dedup, [ [:], [] ] )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam            = UMICOLLAPSE.out.bam             // channel: [ val(meta), path(bam) ]

    bai            = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), path(bai) ]
    csi            = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), path(csi) ]
    dedup_stats    = UMICOLLAPSE.out.log             // channel: [ val(meta), path(stats) ]
    stats          = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat       = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats       = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions       = ch_versions                     // channel: [ path(versions.yml) ]
}
