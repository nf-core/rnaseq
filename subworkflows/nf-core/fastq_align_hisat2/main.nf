include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'
include { Hisat2Aligned           } from './types'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/hisat2/index
    splicesites // channel: /path/to/genome.splicesites.txt
    ch_fasta_fai // channel: [meta, fasta, fai ]
    save_unaligned // val: boolean

    main:
    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN(reads, index, splicesites, save_unaligned)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_orig_bam = HISAT2_ALIGN.out.map { r -> [ r.meta, r.orig_bam ] }
    ch_summary = HISAT2_ALIGN.out.map { r -> [ r.meta, r.hisat2.summary ] }
    ch_fastq = HISAT2_ALIGN.out.filter { r -> r.unmapped }.map { r -> [ r.meta, r.unmapped ] }

    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, ch_fasta_fai)

    ch_results = HISAT2_ALIGN.out
        .map { r -> r + record(aligner: 'hisat2') }
        .join(BAM_SORT_STATS_SAMTOOLS.out.results, by: 'id')

    emit:
    orig_bam = ch_orig_bam // channel: [ val(meta), bam   ]
    summary  = ch_summary // channel: [ val(meta), log   ]
    fastq    = ch_fastq // channel: [ val(meta), fastq ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam // channel: [ val(meta), [ bam ] ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    results  = ch_results // channel: Hisat2Aligned
}
