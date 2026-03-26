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
    reads         // channel: [ val(meta), [ reads ] ]
    index         // channel: /path/to/bowtie2/index/
    fasta_fai     // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // Map reads with Bowtie2
    //
    BOWTIE2_ALIGN(
        reads,
        index.map { index_path -> [ [id: 'genome'], index_path ] },
        [ [:], [] ],    // No fasta needed for BAM output
        params.save_unaligned,  // save_unaligned - enable for downstream analysis of unmapped reads
        false           // sort_bam - we'll sort with samtools for consistency
    )

    ch_orig_bam = BOWTIE2_ALIGN.out.bam
    ch_log = BOWTIE2_ALIGN.out.log

    // Parse alignment rate from log
    ch_percent_mapped = ch_log.map { meta, log_file -> [ meta, getBowtie2PercentMapped(log_file) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta_fai)

    emit:
    orig_bam       = ch_orig_bam                          // channel: [ val(meta), bam ]
    log_final      = ch_log                               // channel: [ val(meta), log ]
    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    index          = BAM_SORT_STATS_SAMTOOLS.out.index    // channel: [ val(meta), [ index ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped                    // channel: [ val(meta), percent_mapped ]
}
