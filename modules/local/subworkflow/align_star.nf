/*
 * Alignment with STAR
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { UNTAR               } from '../process/untar'                                addParams( options: params.index_options    )
include { STAR_GENOMEGENERATE } from '../../nf-core/software/star/genomegenerate/main' addParams( options: params.index_options    )
include { STAR_ALIGN          } from '../../nf-core/software/star/align/main'          addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS   } from '../../nf-core/subworkflow/bam_sort_samtools'     addParams( options: params.samtools_options )

workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/star/index/
    fasta //    file: /path/to/genome.fasta
    gtf   //    file: /path/to/genome.gtf
    
    main:
    /*
     * Uncompress STAR index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = STAR_GENOMEGENERATE ( fasta, gtf ).index
    }

    /*
     * Map reads with STAR
     */
    STAR_ALIGN ( reads, ch_index, gtf )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
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
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
