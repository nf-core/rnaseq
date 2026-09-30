nextflow.enable.types = true

include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'
include { Hisat2Reads; Hisat2Aligned } from './types'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_samples: Channel<Hisat2Reads>
    index: Value<Path>
    splicesites: Value<Path?>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    save_unaligned: Boolean

    main:
    //
    // Map reads with HISAT2
    //
    ch_hisat2 = HISAT2_ALIGN(ch_samples, index, splicesites, save_unaligned)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_sorted = BAM_SORT_STATS_SAMTOOLS(ch_hisat2, ch_fasta, ch_fai)

    ch_results = ch_hisat2
        .map { r -> r + record(aligner: 'hisat2') }
        .join(ch_sorted, by: 'id')

    emit:
    ch_results
}
