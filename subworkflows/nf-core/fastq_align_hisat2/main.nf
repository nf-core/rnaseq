include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/hisat2/index
    splicesites // channel: /path/to/genome.splicesites.txt
    ch_fasta    // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()


    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN ( reads, index, splicesites )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( HISAT2_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)


    emit:
    orig_bam = HISAT2_ALIGN.out.bam                 // channel: [ val(meta), bam   ]
    summary  = HISAT2_ALIGN.out.summary             // channel: [ val(meta), log   ]
    fastq    = HISAT2_ALIGN.out.fastq               // channel: [ val(meta), fastq ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}

