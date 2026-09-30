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

    ch_dedup_bam = UMICOLLAPSE.out
        .filter { r -> r.bam }
        .map { r -> [r.meta, r.bam] }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX(ch_dedup_bam)

    ch_index = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.index] }

    ch_bam_bai_dedup = ch_dedup_bam.join(ch_index, by: [0])

    BAM_STATS_SAMTOOLS(ch_bam_bai_dedup, [[:], [], []])

    ch_results = UMICOLLAPSE.out
        .filter { r -> r.bam }
        .map { r -> [r.id, r] }
        .join(SAMTOOLS_INDEX.out.map { r -> [r.meta.id, r.index] }, by: [0])
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0])
        .map { id, dedup, bai, samtools ->
            record(
                id:          id,
                bam:         dedup.bam,
                bai:         bai,
                dedup_stats: dedup.log,
                samtools:    record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats)
            )
        }

    emit:
    bam         = ch_dedup_bam // channel: [ val(meta), path(bam) ]
    index       = ch_index // channel: [ val(meta), path(index) ]
    dedup_stats = UMICOLLAPSE.out.map { r -> [r.meta, r.log] } // channel: [ val(meta), path(stats) ]
    stats       = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat    = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats    = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    results     = ch_results // channel: UmicollapseDedupBam
}
