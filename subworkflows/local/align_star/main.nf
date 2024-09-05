//
// Alignment with STAR
//

include { STAR_ALIGN              } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN_IGENOMES     } from '../../../modules/local/star_align_igenomes'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'

workflow ALIGN_STAR {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    fasta               // channel: /path/to/fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    ch_orig_bam       = Channel.empty()
    ch_log_final      = Channel.empty()
    ch_log_out        = Channel.empty()
    ch_log_progress   = Channel.empty()
    ch_bam_sorted     = Channel.empty()
    ch_bam_transcript = Channel.empty()
    ch_fastq          = Channel.empty()
    ch_tab            = Channel.empty()
    if (is_aws_igenome) {
        STAR_ALIGN_IGENOMES ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_orig_bam       = STAR_ALIGN_IGENOMES.out.bam
        ch_log_final      = STAR_ALIGN_IGENOMES.out.log_final
        ch_log_out        = STAR_ALIGN_IGENOMES.out.log_out
        ch_log_progress   = STAR_ALIGN_IGENOMES.out.log_progress
        ch_bam_sorted     = STAR_ALIGN_IGENOMES.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN_IGENOMES.out.bam_transcript
        ch_fastq          = STAR_ALIGN_IGENOMES.out.fastq
        ch_tab            = STAR_ALIGN_IGENOMES.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN_IGENOMES.out.versions.first())
    } else {
        STAR_ALIGN ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_orig_bam       = STAR_ALIGN.out.bam
        ch_log_final      = STAR_ALIGN.out.log_final
        ch_log_out        = STAR_ALIGN.out.log_out
        ch_log_progress   = STAR_ALIGN.out.log_progress
        ch_bam_sorted     = STAR_ALIGN.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN.out.bam_transcript
        ch_fastq          = STAR_ALIGN.out.fastq
        ch_tab            = STAR_ALIGN.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( ch_orig_bam, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    orig_bam       = ch_orig_bam                    // channel: [ val(meta), bam            ]
    log_final      = ch_log_final                   // channel: [ val(meta), log_final      ]
    log_out        = ch_log_out                     // channel: [ val(meta), log_out        ]
    log_progress   = ch_log_progress                // channel: [ val(meta), log_progress   ]
    bam_sorted     = ch_bam_sorted                  // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq          = ch_fastq                       // channel: [ val(meta), fastq          ]
    tab            = ch_tab                         // channel: [ val(meta), tab            ]

    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi            = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                    // channel: [ versions.yml ]
}
