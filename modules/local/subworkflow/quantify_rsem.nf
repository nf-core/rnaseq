/*
 * Gene/transcript quantification with RSEM
 */

include { UNTAR                    } from '../process/untar'
include { RSEM_MERGE_COUNTS        } from '../process/rsem_merge_counts'

include { RSEM_PREPAREREFERENCE    } from '../../nf-core/software/rsem/preparereference/main'
include { RSEM_CALCULATEEXPRESSION } from '../../nf-core/software/rsem/calculateexpression/main'
include { BAM_SORT_SAMTOOLS        } from '../../nf-core/subworkflow/bam_sort_samtools'

workflow QUANTIFY_RSEM {
    take:
    reads                       // channel: [ val(meta), [ reads ] ]
    index                       //    file: /path/to/rsem/index/
    fasta                       //    file: /path/to/genome.fasta
    gtf                         //    file: /path/to/genome.gtf
    publish_index_options       //     map: options for publishing index
    preparereference_options    //     map: options for rsem_preparereference module
    calculateexpression_options //     map: options for rsem_calculateexpression module
    samtools_options            //     map: options for bam_sort_samtools subworkflow
    merge_counts_options        //     map: options for merge_counts_rsem module

    main:
    /*
     * Uncompress RSEM index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, publish_index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = RSEM_PREPAREREFERENCE ( fasta, gtf, preparereference_options ).index
    }

    /*
     * Quantify reads with RSEM
     */
    RSEM_CALCULATEEXPRESSION ( reads, ch_index, calculateexpression_options )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( RSEM_CALCULATEEXPRESSION.out.bam_star, samtools_options )

    /*
     * Merge counts across samples
     */
    RSEM_MERGE_COUNTS (
        RSEM_CALCULATEEXPRESSION.out.counts_gene.collect{it[1]},
        RSEM_CALCULATEEXPRESSION.out.counts_transcript.collect{it[1]},
        merge_counts_options
    )

    emit:
    counts_gene              = RSEM_CALCULATEEXPRESSION.out.counts_gene       // channel: [ val(meta), counts ]
    counts_transcript        = RSEM_CALCULATEEXPRESSION.out.counts_transcript // channel: [ val(meta), counts ]
    stat                     = RSEM_CALCULATEEXPRESSION.out.stat              // channel: [ val(meta), stat ]
    logs                     = RSEM_CALCULATEEXPRESSION.out.logs              // channel: [ val(meta), logs ]
    bam_star                 = RSEM_CALCULATEEXPRESSION.out.bam_star          // channel: [ val(meta), bam ]
    bam_genome               = RSEM_CALCULATEEXPRESSION.out.bam_genome        // channel: [ val(meta), bam ]
    bam_transcript           = RSEM_CALCULATEEXPRESSION.out.bam_transcript    // channel: [ val(meta), bam ]
    rsem_version             = RSEM_CALCULATEEXPRESSION.out.version           //    path: *.version.txt

    bam                      = BAM_SORT_SAMTOOLS.out.bam                      // channel: [ val(meta), [ bam ] ]
    bai                      = BAM_SORT_SAMTOOLS.out.bai                      // channel: [ val(meta), [ bai ] ]
    stats                    = BAM_SORT_SAMTOOLS.out.stats                    // channel: [ val(meta), [ stats ] ]
    flagstat                 = BAM_SORT_SAMTOOLS.out.flagstat                 // channel: [ val(meta), [ flagstat ] ]
    idxstats                 = BAM_SORT_SAMTOOLS.out.idxstats                 // channel: [ val(meta), [ idxstats ] ]
    samtools_version         = BAM_SORT_SAMTOOLS.out.samtools_version         //    path: *.version.txt

    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.tsv
    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.tsv
    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.tsv
    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.tsv
}
