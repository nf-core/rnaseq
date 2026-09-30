nextflow.enable.types = true

//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { RawBams; Bam } from '../../../modules/nf-core/types'
include { SamtoolsIndexResult } from '../../../modules/nf-core/samtools/index/main'
include { SamtoolsSortResult } from '../../../modules/nf-core/samtools/sort/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam: Channel<RawBams>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    def ch_sorted: Channel<SamtoolsSortResult> = SAMTOOLS_SORT(ch_bam, ch_fasta, ch_fai, '')
        .filter { r -> r.bam != null }

    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_sorted)
    ch_sorted_indexed = ch_sorted.join(ch_index, by: 'id', remainder: true)
    ch_sorted_indexed.subscribe { r ->
        if( r.bam == null || r.bai == null ) {
            error "Sample '${r.id}' is missing its samtools index result"
        }
    }
    ch_indexed = ch_sorted_indexed.filter { r -> r.bam != null && r.bai != null }

    // SAMTOOLS_SORT also carries cram, sam, csi and crai fields; dropping them keeps them from
    // overwriting same-named fields when a caller joins this result onto its own record.
    ch_results = ch_indexed
        .join(BAM_STATS_SAMTOOLS(ch_indexed, ch_fasta, ch_fai), by: 'id')
        .map { r -> record(id: r.id, meta: r.meta, bam: r.bam, bai: r.bai, samtools: r.samtools) }

    emit:
    ch_results
}
