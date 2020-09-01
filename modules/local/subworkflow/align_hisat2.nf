/*
 * Alignment with HISAT2
 */

include { UNTAR                     } from '../process/untar'
include { HISAT2_EXTRACTSPLICESITES } from '../process/hisat2_extractsplicesites'
include { HISAT2_BUILD              } from '../process/hisat2_build'
include { HISAT2_ALIGN              } from '../process/hisat2_align'
include { BAM_SORT_SAMTOOLS         } from '../../nf-core/subworkflow/bam_sort_samtools'

workflow ALIGN_HISAT2 {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    index            //    file: /path/to/star/index/
    fasta            //    file: /path/to/genome.fasta
    gtf              //    file: /path/to/genome.gtf
    splicesites      //    file: /path/to/genome.splicesites.txt
    index_options    //     map: options for hisat2_build module
    align_options    //     map: options for hisat2_align module
    samtools_options //     map: options for bam_sort_samtools subworkflow

    main:
    /*
     * Uncompress HISAT2 index or generate from scratch if required
    */
    if (!splicesites) {
        ch_splicesites = HISAT2_EXTRACTSPLICESITES ( gtf, index_options ).txt
    } else {
        ch_splicesites = file(splicesites)
    }
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = HISAT2_BUILD ( fasta, gtf, ch_splicesites, index_options ).index
    }

    /*
     * Map reads with HISAT2
     */
    HISAT2_ALIGN ( reads, ch_index, ch_splicesites, align_options )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( HISAT2_ALIGN.out.bam, samtools_options )

    emit:
    orig_bam         = HISAT2_ALIGN.out.bam     // channel: [ val(meta), bam   ]
    summary          = HISAT2_ALIGN.out.summary // channel: [ val(meta), log   ]
    fastq            = HISAT2_ALIGN.out.fastq   // channel: [ val(meta), fastq ]
    hisat2_version   = HISAT2_ALIGN.out.version //    path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam              // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai              // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats            // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat         // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats         // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.samtools_version //    path: *.version.txt
}
