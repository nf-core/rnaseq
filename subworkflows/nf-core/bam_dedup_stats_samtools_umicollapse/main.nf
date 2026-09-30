nextflow.enable.types = true

//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { Bam             } from '../../local/types'
include { UmicollapseDedupBam } from './types'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam_bai: Channel<Bam>

    main:
    ch_no_fasta = channel.value(tuple([:], null, null))

    //
    // umicollapse in bam mode (thus hardcode mode input to 'bam')
    //
    ch_dedup = UMICOLLAPSE(ch_bam_bai, 'bam')
        .filter { r -> r.bam != null }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_indexed = ch_dedup.join(SAMTOOLS_INDEX(ch_dedup), by: 'id')

    ch_results = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, ch_no_fasta), by: 'id')

    emit:
    ch_results
}
