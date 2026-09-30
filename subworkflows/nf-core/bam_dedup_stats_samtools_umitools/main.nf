nextflow.enable.types = true

//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'
include { Bam; UmitoolsDedupBam } from '../../../modules/nf-core/types'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_bai: Channel<Bam>
    val_get_dedup_stats: Boolean
    val_primary_only: Boolean

    main:
    //
    // Optionally filter to primary alignments before deduplication
    //
    ch_dedup_input = ch_bam_bai
    if (val_primary_only) {
        ch_primary = SAMTOOLS_VIEW_PRIMARY(
            ch_bam_bai,
            null,
            null,
            null,
            null,
            '',
        ).filter { r -> r.bam != null }

        ch_dedup_input = ch_primary.join(SAMTOOLS_INDEX_PRIMARY(ch_primary), by: 'id')
    }

    //
    // UMI-tools dedup
    //
    ch_dedup = UMITOOLS_DEDUP(ch_dedup_input, val_get_dedup_stats)

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_indexed = ch_dedup.join(SAMTOOLS_INDEX(ch_dedup), by: 'id')

    ch_results = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, null, null), by: 'id')

    emit:
    ch_results
}
