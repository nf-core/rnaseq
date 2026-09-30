//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'
include { UmitoolsDedupBam                         } from './types'

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
            [[:], [], []],
            [[:], []],
            [[:], []],
            '',
        )

        ch_primary_bam = SAMTOOLS_VIEW_PRIMARY.out
            .filter { r -> r.bam }
            .map { r -> [r.meta, r.bam] }

        SAMTOOLS_INDEX_PRIMARY(ch_primary_bam)

        ch_dedup_input = ch_primary_bam.join(SAMTOOLS_INDEX_PRIMARY.out.map { r -> [r.meta, r.index] }, by: [0])
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
    ch_dedup_bam = UMITOOLS_DEDUP.out.map { r -> [r.meta, r.bam] }

    SAMTOOLS_INDEX(ch_dedup_bam)

    ch_index = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.index] }

    ch_bam_bai_dedup = ch_dedup_bam.join(ch_index, by: [0])

    BAM_STATS_SAMTOOLS(ch_bam_bai_dedup, [[:], [], []])

    // umi_tools writes all three stats tables together, so `tsv` is either
    // complete or absent.
    ch_results = UMITOOLS_DEDUP.out
        .map { r -> [r.id, r] }
        .join(SAMTOOLS_INDEX.out.map { r -> [r.meta.id, r.index] }, by: [0])
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0])
        .map { id, dedup, bai, samtools ->
            record(
                id:        id,
                bam:       dedup.bam,
                bai:       bai,
                dedup_log: dedup.log,
                samtools:  record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats),
                tsv:       [dedup.tsv_edit_distance, dedup.tsv_per_umi, dedup.tsv_umi_per_position].any()
                    ? record(edit_distance: dedup.tsv_edit_distance, per_umi: dedup.tsv_per_umi, umi_per_position: dedup.tsv_umi_per_position)
                    : null
            )
        }

    emit:
    bam                  = ch_dedup_bam // channel: [ val(meta), path(bam) ]
    deduplog             = UMITOOLS_DEDUP.out.map { r -> [r.meta, r.log] } // channel: [ val(meta), path(log) ]
    tsv_edit_distance    = UMITOOLS_DEDUP.out.filter { r -> r.tsv_edit_distance }.map { r -> [r.meta, r.tsv_edit_distance] } // channel: [ val(meta), path(tsv) ]
    tsv_per_umi          = UMITOOLS_DEDUP.out.filter { r -> r.tsv_per_umi }.map { r -> [r.meta, r.tsv_per_umi] } // channel: [ val(meta), path(tsv) ]
    tsv_umi_per_position = UMITOOLS_DEDUP.out.filter { r -> r.tsv_umi_per_position }.map { r -> [r.meta, r.tsv_umi_per_position] } // channel: [ val(meta), path(tsv) ]
    index                = ch_index // channel: [ val(meta), path(index) ]
    stats                = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat             = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats             = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    results              = ch_results // channel: UmitoolsDedupBam
}
