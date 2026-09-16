//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai/csi) ]
    val_get_dedup_stats // boolean: true/false
    val_primary_only // boolean: true/false

    main:

    //
    // Optionally filter to primary alignments before deduplication
    //
    if (val_primary_only) {
        SAMTOOLS_VIEW_PRIMARY(
            ch_bam_bai,
            [[], [], []],
            [[], []],
            [[], []],
            [],
        )

        SAMTOOLS_INDEX_PRIMARY(SAMTOOLS_VIEW_PRIMARY.out.bam)

        ch_dedup_input = SAMTOOLS_VIEW_PRIMARY.out.bam.join(SAMTOOLS_INDEX_PRIMARY.out.index, by: [0])
    }
    else {
        ch_dedup_input = ch_bam_bai
    }

    //
    // UMI-tools dedup
    //
    UMITOOLS_DEDUP(ch_dedup_input, val_get_dedup_stats)

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX(UMITOOLS_DEDUP.out.bam)

    ch_bam_bai_dedup = UMITOOLS_DEDUP.out.bam.join(SAMTOOLS_INDEX.out.index, by: [0])

    BAM_STATS_SAMTOOLS(ch_bam_bai_dedup, [[:], [], []])

    emit:
    bam                  = UMITOOLS_DEDUP.out.bam // channel: [ val(meta), path(bam) ]
    deduplog             = UMITOOLS_DEDUP.out.log // channel: [ val(meta), path(log) ]
    tsv_edit_distance    = UMITOOLS_DEDUP.out.tsv_edit_distance // channel: [ val(meta), path(tsv) ]
    tsv_per_umi          = UMITOOLS_DEDUP.out.tsv_per_umi // channel: [ val(meta), path(tsv) ]
    tsv_umi_per_position = UMITOOLS_DEDUP.out.tsv_umi_per_position // channel: [ val(meta), path(tsv) ]
    index                = SAMTOOLS_INDEX.out.index // channel: [ val(meta), path(index) ]
    stats                = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat             = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats             = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
}
