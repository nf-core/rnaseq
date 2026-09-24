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

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX(UMICOLLAPSE.out.bam)

    ch_bam_bai_dedup = UMICOLLAPSE.out.bam.join(SAMTOOLS_INDEX.out.index, by: [0])

    BAM_STATS_SAMTOOLS(ch_bam_bai_dedup, [[:], [], []])

    ch_results = ch_bam_bai_dedup
        .join(UMICOLLAPSE.out.log, by: [0])
        .map { meta, bam, bai, dedup_stats -> [meta.id, bam, bai, dedup_stats] }
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0])
        .map { id, bam, bai, dedup_stats, samtools ->
            record(
                id:          id,
                bam:         bam,
                bai:         bai,
                dedup_stats: dedup_stats,
                samtools:    record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats)
            )
        }

    emit:
    bam         = UMICOLLAPSE.out.bam // channel: [ val(meta), path(bam) ]
    index       = SAMTOOLS_INDEX.out.index // channel: [ val(meta), path(index) ]
    dedup_stats = UMICOLLAPSE.out.log // channel: [ val(meta), path(stats) ]
    stats       = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat    = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats    = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    results     = ch_results // channel: UmicollapseDedupBam
}
