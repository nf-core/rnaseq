nextflow.enable.types = true

//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { SamtoolsIndexResult; UmicollapseResult; Bam; UmicollapseDedupBam } from '../../../modules/nf-core/types'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam_bai: Channel<Bam>

    main:
    //
    // umicollapse in bam mode (thus hardcode mode input to 'bam')
    //
    def ch_dedup: Channel<UmicollapseResult> = UMICOLLAPSE(ch_bam_bai, 'bam')
        .filter { r -> r.bam != null }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_dedup)
    ch_indexed = ch_dedup.join(ch_index, by: 'id')

    def ch_results: Channel<UmicollapseDedupBam> = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, null, null), by: 'id')

    emit:
    ch_results
}
