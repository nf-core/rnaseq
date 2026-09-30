//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { UmicollapseDedupBam } from './types'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai/csi) ]

    main:
    //
    // umicollapse in bam mode (thus hardcode mode input channel to 'bam')
    //
    UMICOLLAPSE(ch_bam_bai, channel.value('bam'))

    ch_dedup = UMICOLLAPSE.out.filter { r -> r.bam }
    ch_dedup_bam = ch_dedup.map { r -> [r.meta, r.bam] }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX(ch_dedup)

    ch_index = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.bai] }

    ch_bam_bai_dedup = ch_dedup.join(SAMTOOLS_INDEX.out, by: 'id')

    BAM_STATS_SAMTOOLS(ch_bam_bai_dedup, [[:], [], []])

    ch_results = ch_bam_bai_dedup
        .join(BAM_STATS_SAMTOOLS.out, by: 'id')
        .map { r ->
            record(
                id:          r.id,
                bam:         r.bam,
                bai:         r.bai,
                dedup_stats: r.log,
                samtools:    r.samtools
            )
        }

    emit:
    bam         = ch_dedup_bam // channel: [ val(meta), path(bam) ]
    index       = ch_index // channel: [ val(meta), path(index) ]
    dedup_stats = UMICOLLAPSE.out.map { r -> [r.meta, r.log] } // channel: [ val(meta), path(stats) ]
    stats       = BAM_STATS_SAMTOOLS.out.map { r -> [r.meta, r.samtools.stats] } // channel: [ val(meta), path(stats) ]
    flagstat    = BAM_STATS_SAMTOOLS.out.map { r -> [r.meta, r.samtools.flagstat] } // channel: [ val(meta), path(flagstat) ]
    idxstats    = BAM_STATS_SAMTOOLS.out.map { r -> [r.meta, r.samtools.idxstats] } // channel: [ val(meta), path(idxstats) ]
    results     = ch_results // channel: UmicollapseDedupBam
}
