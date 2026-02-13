include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/hisat2/index
    splicesites // channel: /path/to/genome.splicesites.txt
    ch_fasta    // channel: [ fasta ]

    main:

    ch_versions = channel.empty()


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

    // Extract individual fields from SamtoolsResult record
    ch_samtools = BAM_SORT_STATS_SAMTOOLS.out.result

    emit:
    orig_bam = HISAT2_ALIGN.out.bam                      // channel: [ val(meta), bam   ]
    summary  = HISAT2_ALIGN.out.summary                   // channel: [ val(meta), log   ]
    fastq    = HISAT2_ALIGN.out.fastq                     // channel: [ val(meta), fastq ]

    bam      = ch_samtools.map { r -> [r.meta, r.bam] }      // channel: [ val(meta), [ bam ] ]
    bai      = ch_samtools.map { r -> [r.meta, r.bai] }      // channel: [ val(meta), [ bai ] ]
    csi      = ch_samtools.map { r -> [r.meta, r.csi] }      // channel: [ val(meta), [ csi ] ]
    stats    = ch_samtools.map { r -> [r.meta, r.stats] }    // channel: [ val(meta), [ stats ] ]
    flagstat = ch_samtools.map { r -> [r.meta, r.flagstat] }  // channel: [ val(meta), [ flagstat ] ]
    idxstats = ch_samtools.map { r -> [r.meta, r.idxstats] }  // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                                 // channel: [ versions.yml ]
}
