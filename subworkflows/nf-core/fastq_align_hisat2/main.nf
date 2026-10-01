nextflow.enable.types = true

include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'
include { ReadsInput; Bam; Hisat2Aligned } from '../../../modules/nf-core/types'
include { Hisat2AlignResult } from '../../../modules/nf-core/hisat2/align/main'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_samples: Channel<ReadsInput>
    index: Value<Path>
    splicesites: Value<Path?>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    save_unaligned: Boolean
    index_args: String // samtools index options

    main:
    //
    // Map reads with HISAT2
    //
    def ch_hisat2: Channel<Hisat2AlignResult> = HISAT2_ALIGN(ch_samples, index, splicesites, save_unaligned)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_sorted: Channel<Bam> = BAM_SORT_STATS_SAMTOOLS(ch_hisat2, ch_fasta, ch_fai, index_args)

    ch_results = ch_hisat2
        .map { r -> r + record(aligner: 'hisat2') }
        .join(ch_sorted, by: 'id')

    emit:
    ch_results
}
