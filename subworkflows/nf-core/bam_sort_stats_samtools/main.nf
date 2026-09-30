nextflow.enable.types = true

//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { Bam; RawBams } from '../../../modules/nf-core/types'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam: Channel<RawBams>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    ch_sorted = SAMTOOLS_SORT(ch_bam, ch_fasta, ch_fai, '')
        .filter { r -> r.bam != null }

    ch_indexed = ch_sorted.join(SAMTOOLS_INDEX(ch_sorted), by: 'id')

    // SAMTOOLS_SORT also carries cram, sam, csi and crai fields; dropping them keeps them from
    // overwriting same-named fields when a caller joins this result onto its own record.
    ch_results = ch_indexed
        .join(BAM_STATS_SAMTOOLS(ch_indexed, ch_fasta, ch_fai), by: 'id')
        .map { r -> record(id: r.id, meta: r.meta, bam: r.bam, bai: r.bai, samtools: r.samtools) }

    emit:
    ch_results
}
