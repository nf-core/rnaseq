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
include { Reads; Bowtie2Aligned   } from './types'

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
    ch_samples: Channel<Reads>
    index: Value<Path> // /path/to/bowtie2/index/
    fasta_fai: Value<Tuple<Map, Path?, Path?>>

    main:

    //
    // Map reads with Bowtie2
    //
    ch_bowtie2 = BOWTIE2_ALIGN(
        ch_samples,
        index.map { index_path -> tuple([id: 'genome'], index_path) },
        tuple([:], null),       // No fasta needed for BAM output
        params.save_unaligned,  // save_unaligned - enable for downstream analysis of unmapped reads
        false                   // sort_bam - we'll sort with samtools for consistency
    ).filter { r -> !r.orig_bam.isEmpty() }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_sorted = BAM_SORT_STATS_SAMTOOLS(
        ch_bowtie2.map { r -> record(id: r.id, meta: r.meta, bam: r.orig_bam) },
        fasta_fai
    )

    // Parse alignment rate from log
    ch_results = ch_bowtie2
        .map { r -> r + record(aligner: 'bowtie2', percent_mapped: getBowtie2PercentMapped(r.bowtie2.log)) }
        .join(ch_sorted, by: 'id')

    emit:
    ch_results // channel: Bowtie2Aligned
}
