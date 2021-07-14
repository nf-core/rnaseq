//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

params.markduplicates_options = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main' addParams( options: params.markduplicates_options )
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'        addParams( options: params.samtools_index_options )
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'                                     addParams( options: params.samtools_stats_options )

workflow MARK_DUPLICATES_PICARD {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    //
    // Picard MarkDuplicates
    //
    PICARD_MARKDUPLICATES ( bam )

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )

    PICARD_MARKDUPLICATES.out.bam
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
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai )

    emit:
    bam              = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), [ bam ] ]
    metrics          = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), [ metrics ] ]
    picard_version   = PICARD_MARKDUPLICATES.out.version //    path: *.version.txt

    bai              = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    csi              = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_INDEX.out.version        //    path: *.version.txt
}
