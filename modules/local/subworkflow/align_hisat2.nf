/*
 * Alignment with HISAT2
 */

params.align_options    = [:]
params.samtools_options = [:]

include { HISAT2_ALIGN      } from '../../nf-core/software/hisat2/align/main'    addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from '../../nf-core/subworkflow/bam_sort_samtools' addParams( options: params.samtools_options )

workflow ALIGN_HISAT2 {
    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/star/index/
    splicesites // channel: /path/to/genome.splicesites.txt
    
    main:
    /*
     * Map reads with HISAT2
     */
    HISAT2_ALIGN ( reads, index, splicesites )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( HISAT2_ALIGN.out.bam )

    emit:
    orig_bam         = HISAT2_ALIGN.out.bam           // channel: [ val(meta), bam   ]
    summary          = HISAT2_ALIGN.out.summary       // channel: [ val(meta), log   ]
    fastq            = HISAT2_ALIGN.out.fastq         // channel: [ val(meta), fastq ]
    hisat2_version   = HISAT2_ALIGN.out.version       //    path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
