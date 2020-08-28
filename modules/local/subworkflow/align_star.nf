/*
 * Alignment with STAR
 */

include { UNTAR               } from '../process/untar'
include { STAR_GENOMEGENERATE } from '../process/star_genomegenerate'
include { STAR_ALIGN          } from '../process/star_align'

workflow ALIGN_STAR {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index         //    file: /path/to/star/index/
    fasta         //    file: /path/to/genome.fasta
    gtf           //    file: /path/to/genome.gtf
    index_options //     map: options for star_genomegenerate module
    align_options //     map: options for star_align module

    main:
    /*
     * Uncompress STAR index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = STAR_GENOMEGENERATE ( fasta, gtf, index_options ).index
    }

    /*
     * Map reads with STAR
     */
    STAR_ALIGN ( reads, ch_index, gtf, align_options )

    emit:
    bam            = STAR_ALIGN.out.bam            // channel: [ val(meta), bam            ]
    log_final      = STAR_ALIGN.out.log_final      // channel: [ val(meta), log_final      ]
    log_out        = STAR_ALIGN.out.log_out        // channel: [ val(meta), log_out        ]
    log_progress   = STAR_ALIGN.out.log_progress   // channel: [ val(meta), log_progress   ]
    bam_sorted     = STAR_ALIGN.out.bam_sorted     // channel: [ val(meta), bam_sorted     ]
    bam_transcript = STAR_ALIGN.out.bam_transcript // channel: [ val(meta), bam_transcript ]
    fastq          = STAR_ALIGN.out.fastq          // channel: [ val(meta), fastq          ]
    tab            = STAR_ALIGN.out.tab            // channel: [ val(meta), tab            ]
    version        = STAR_ALIGN.out.version        // path: *.version.txt
}
