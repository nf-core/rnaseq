//
// Alignment with STAR
//

params.align_options          = [:]
params.samtools_sort_options  = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { STAR_ALIGN        } from '../../modules/nf-core/software/star/align/main' addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'                            addParams( sort_options: params.samtools_sort_options, index_options: params.samtools_index_options, stats_options: params.samtools_stats_options )

workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    gtf   // channel: /path/to/genome.gtf

    main:

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( reads, index, gtf )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( STAR_ALIGN.out.bam )

    emit:
    orig_bam         = STAR_ALIGN.out.bam             // channel: [ val(meta), bam            ]
    log_final        = STAR_ALIGN.out.log_final       // channel: [ val(meta), log_final      ]
    log_out          = STAR_ALIGN.out.log_out         // channel: [ val(meta), log_out        ]
    log_progress     = STAR_ALIGN.out.log_progress    // channel: [ val(meta), log_progress   ]
    bam_sorted       = STAR_ALIGN.out.bam_sorted      // channel: [ val(meta), bam_sorted     ]
    bam_transcript   = STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), bam_transcript ]
    fastq            = STAR_ALIGN.out.fastq           // channel: [ val(meta), fastq          ]
    tab              = STAR_ALIGN.out.tab             // channel: [ val(meta), tab            ]
    star_version     = STAR_ALIGN.out.version         // path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi              = BAM_SORT_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
