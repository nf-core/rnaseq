nextflow.enable.types = true

//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'
include { BamBaiInput; UmitoolsDedupBam } from '../../../modules/nf-core/types'
include { SamtoolsIndexResult } from '../../../modules/nf-core/samtools/index/main'
include { SamtoolsViewResult } from '../../../modules/nf-core/samtools/view/main'
include { UmitoolsDedupResult } from '../../../modules/nf-core/umitools/dedup/main'

//
// UMI-tools dedup options
//
def umitoolsDedupArgs(grouping_method, umi_separator, meta) {
    return [
        meta.single_end ? '' : '--unpaired-reads=discard --chimeric-pairs=discard',
        grouping_method ? "--method='${grouping_method}'" : '',
        umi_separator   ? "--umi-separator='${umi_separator}'" : ''
    ].join(' ').trim()
}

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_bai: Channel<BamBaiInput>
    val_get_dedup_stats: Boolean
    val_primary_only: Boolean
    index_args: String // samtools index options
    grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:

    //
    // Optionally filter to primary alignments before deduplication
    //
    if (val_primary_only) {
        def ch_primary: Channel<SamtoolsViewResult> = SAMTOOLS_VIEW_PRIMARY(
            ch_bam_bai,
            null,
            null,
            null,
            null,
            '',
        ).filter { r -> r.bam != null }

        ch_dedup_input = ch_primary.join(SAMTOOLS_INDEX_PRIMARY(ch_primary.map { r -> r + record(args: index_args) }), by: 'id')
    }
    else {
        ch_dedup_input = ch_bam_bai
    }

    //
    // UMI-tools dedup
    //
    def ch_dedup: Channel<UmitoolsDedupResult> = UMITOOLS_DEDUP(
        ch_dedup_input.map { r -> r + record(args: umitoolsDedupArgs(grouping_method, umi_separator, r.meta)) },
        val_get_dedup_stats
    )

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_dedup.map { r -> r + record(args: index_args) })
    ch_indexed = ch_dedup.join(ch_index, by: 'id')

    ch_results = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, null, null), by: 'id')

    emit:
    ch_results
}
