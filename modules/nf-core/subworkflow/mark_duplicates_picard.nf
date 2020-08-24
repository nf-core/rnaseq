/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES } from '../software/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../software/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'

workflow MARK_DUPLICATES_PICARD {
    take:
    ch_bam                 // channel: [ val(meta), [ bam ] ]
    markduplicates_options //     map: options for picard MarkDuplicates module
    samtools_options       //     map: options for SAMTools modules

    main:
    PICARD_MARKDUPLICATES(ch_bam, markduplicates_options)
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam, samtools_options)
    BAM_STATS_SAMTOOLS(PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]), samtools_options)

    emit:
    bam = PICARD_MARKDUPLICATES.out.bam                // channel: [ val(meta), [ bam ] ]
    metrics = PICARD_MARKDUPLICATES.out.metrics        // channel: [ val(meta), [ metrics ] ]
    bai = SAMTOOLS_INDEX.out.bai                       // channel: [ val(meta), [ bai ] ]
    stats = BAM_STATS_SAMTOOLS.out.stats               // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat         // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats         // channel: [ val(meta), [ idxstats ] ]
    picard_version = PICARD_MARKDUPLICATES.out.version //    path: *.version.txt
    samtools_version = SAMTOOLS_INDEX.out.version      //    path: *.version.txt
}
