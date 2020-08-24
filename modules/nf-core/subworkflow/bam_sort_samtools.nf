/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { SAMTOOLS_SORT      } from '../software/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../software/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools'

workflow BAM_SORT_SAMTOOLS {
    take:
    ch_bam           // channel: [ val(meta), [ bam ] ]
    samtools_options //     map: options for SAMTools modules

    main:
    SAMTOOLS_SORT(ch_bam, samtools_options)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam, samtools_options)
    BAM_STATS_SAMTOOLS(SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), samtools_options)

    emit:
    bam = SAMTOOLS_SORT.out.bam                  // channel: [ val(meta), [ bam ] ]
    bai = SAMTOOLS_INDEX.out.bai                 // channel: [ val(meta), [ bai ] ]
    stats = BAM_STATS_SAMTOOLS.out.stats         // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_SORT.out.version //    path: *.version.txt
}
