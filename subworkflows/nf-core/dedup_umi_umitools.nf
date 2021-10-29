//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

// params.umitools_options       = [:]
// params.samtools_index_options = [:]
// params.samtools_stats_options = [:]

include { UMITOOLS_DEDUP     } from '../../modules/nf-core/modules/umitools/dedup/main' // addParams( options: params.umitools_options       )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' // addParams( options: params.samtools_index_options )
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools'                              // addParams( options: params.samtools_stats_options )

workflow DEDUP_UMI_UMITOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai/csi ] ]

    main:

    ch_versions = Channel.empty()

    //
    // UMI-tools dedup
    //
    UMITOOLS_DEDUP ( bam_bai )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    UMITOOLS_DEDUP.out.bam
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
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam              = UMITOOLS_DEDUP.out.bam          // channel: [ val(meta), [ bam ] ]

    bai              = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi              = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                     // channel: [ versions.yml ]
}
