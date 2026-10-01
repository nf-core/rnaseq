nextflow.enable.types = true

//
// Alignment with Bowtie2
//
// Bowtie2 is a splice-unaware aligner, suitable for prokaryotic RNA-seq data.
// It aligns reads to the transcriptome (not genome) for Salmon quantification.
// The Bowtie2 index is built from transcript FASTA, enabling alignment-based
// Salmon quantification similar to the STAR transcriptome BAM workflow.
//

include { BOWTIE2_ALIGN           } from '../../../modules/nf-core/bowtie2/align'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'
include { ReadsInput; Bowtie2Aligned; Bam } from '../../../modules/nf-core/types'
include { Bowtie2AlignResult } from '../../../modules/nf-core/bowtie2/align/main'

//
// Function that parses and returns the alignment rate from the Bowtie2 log output
//
def getBowtie2PercentMapped(align_log) {
    def percent_aligned = 0
    def pattern = /(\d+\.\d+)% overall alignment rate/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }
    return percent_aligned
}

workflow ALIGN_BOWTIE2 {
    take:
    ch_samples: Channel<ReadsInput>
    index: Value<Path> // /path/to/bowtie2/index/
    fasta: Value<Path?>
    fai: Value<Path?>
    save_unaligned: Boolean
    index_args: String // samtools index options

    main:

    //
    // Map reads with Bowtie2
    //
    def ch_bowtie2: Channel<Bowtie2AlignResult> = BOWTIE2_ALIGN(
        ch_samples,
        index,
        null,                   // no fasta needed for BAM output
        save_unaligned,
        false                   // sort_bam - we'll sort with samtools for consistency
    ).filter { r -> !r.raw_bams.isEmpty() }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_sorted: Channel<Bam> = BAM_SORT_STATS_SAMTOOLS(ch_bowtie2, fasta, fai, index_args)

    // The BAM is aligned to the transcriptome, and its unsorted form is what Salmon quantifies:
    // a coordinate-sorted BAM breaks paired-end quantification.
    ch_bowtie2_sorted = ch_bowtie2
        .map { r ->
            r + record(
                aligner:           'bowtie2',
                percent_mapped:    getBowtie2PercentMapped(r.bowtie2.log),
                transcriptome_bam: r.raw_bams[0]
            )
        }
        .join(ch_sorted, by: 'id', remainder: true)
    ch_bowtie2_sorted.subscribe { r ->
        if( r.bowtie2 == null || r.samtools == null ) {
            error "Sample '${r.id}' is missing its Bowtie2 sorted BAM result"
        }
    }
    ch_results = ch_bowtie2_sorted.filter { r -> r.bowtie2 != null && r.samtools != null }

    emit:
    ch_results
}
