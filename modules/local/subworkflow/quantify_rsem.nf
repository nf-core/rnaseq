/*
 * Gene/transcript quantification with RSEM
 */

params.index_options               = [:]
params.preparereference_options    = [:]
params.calculateexpression_options = [:]
params.samtools_options            = [:]
params.merge_counts_options        = [:]

include { UNTAR                    } from '../process/untar'                                     addParams( options: params.index_options               )
include { RSEM_PREPAREREFERENCE    } from '../../nf-core/software/rsem/preparereference/main'    addParams( options: params.preparereference_options    )
include { RSEM_CALCULATEEXPRESSION } from '../../nf-core/software/rsem/calculateexpression/main' addParams( options: params.calculateexpression_options )
include { BAM_SORT_SAMTOOLS        } from '../../nf-core/subworkflow/bam_sort_samtools'          addParams( options: params.samtools_options            )
include { RSEM_MERGE_COUNTS        } from '../process/rsem_merge_counts'                         addParams( options: params.merge_counts_options        )

workflow QUANTIFY_RSEM {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/rsem/index/
    fasta //    file: /path/to/genome.fasta
    gtf   //    file: /path/to/genome.gtf

    main:
    /*
     * Uncompress RSEM index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = RSEM_PREPAREREFERENCE ( fasta, gtf ).index
    }

    /*
     * Quantify reads with RSEM
     */
    RSEM_CALCULATEEXPRESSION ( reads, ch_index )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( RSEM_CALCULATEEXPRESSION.out.bam_star )

    /*
     * Merge counts across samples
     */
    RSEM_MERGE_COUNTS (
        RSEM_CALCULATEEXPRESSION.out.counts_gene.collect{it[1]},       // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        RSEM_CALCULATEEXPRESSION.out.counts_transcript.collect{it[1]}
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
    samtools_version         = BAM_SORT_SAMTOOLS.out.version                  //    path: *.version.txt

    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.tsv
    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.tsv
    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.tsv
    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.tsv
}
