//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record UmitoolsDedupStats {
    edit_distance:    Path
    per_umi:          Path
    umi_per_position: Path
}

record UmitoolsDedupBam {
    id:        String
    bam:       Path
    bai:       Path
    dedup_log: Path
    samtools:  SamtoolsStatsFiles
    tsv:       UmitoolsDedupStats?
}

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

    // umi_tools writes all three stats tables together, so `tsv` is either
    // complete or absent.
    ch_results = ch_bam_bai_dedup
        .join(UMITOOLS_DEDUP.out.log, by: [0])
        .join(UMITOOLS_DEDUP.out.tsv_edit_distance, by: [0], remainder: true)
        .join(UMITOOLS_DEDUP.out.tsv_per_umi, by: [0], remainder: true)
        .join(UMITOOLS_DEDUP.out.tsv_umi_per_position, by: [0], remainder: true)
        .map { meta, bam, bai, dedup_log, edit_distance, per_umi, umi_per_position ->
            [meta.id, bam, bai, dedup_log, edit_distance, per_umi, umi_per_position]
        }
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0])
        .map { id, bam, bai, dedup_log, edit_distance, per_umi, umi_per_position, samtools ->
            record(
                id:        id,
                bam:       bam,
                bai:       bai,
                dedup_log: dedup_log,
                samtools:  record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats),
                tsv:       [edit_distance, per_umi, umi_per_position].any()
                    ? record(edit_distance: edit_distance, per_umi: per_umi, umi_per_position: umi_per_position)
                    : null
            )
        }

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
    results              = ch_results // channel: UmitoolsDedupBam
}
