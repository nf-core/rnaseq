/*
 * Gene/transcript quantification with RSEM
 */

include { UNTAR                    } from '../process/untar'
include { RSEM_MERGE_COUNTS        } from '../process/rsem_merge_counts'
include { RSEM_PREPAREREFERENCE    } from '../../nf-core/software/rsem/preparereference/main'
include { RSEM_CALCULATEEXPRESSION } from '../../nf-core/software/rsem/calculateexpression/main'

workflow QUANTIFY_RSEM {
    take:
    reads                // channel: [ val(meta), [ reads ] ]
    index                //    file: /path/to/rsem/index/
    fasta                //    file: /path/to/genome.fasta
    gtf                  //    file: /path/to/genome.gtf
    index_options        //     map: options for rsem_preparereference module
    align_options        //     map: options for rsem_calculateexpression module
    merge_counts_options //     map: options for merge_counts_rsem module

    main:
    /*
     * Uncompress RSEM index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = RSEM_PREPAREREFERENCE ( fasta, gtf, index_options ).index
    }

    /*
     * Quantify reads with RSEM
     */
    RSEM_CALCULATEEXPRESSION ( reads, ch_index, align_options )

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
    bam_star                 = RSEM_CALCULATEEXPRESSION.out.bam_star          // channel: [ val(meta), bam ]
    bam_genome               = RSEM_CALCULATEEXPRESSION.out.bam_genome        // channel: [ val(meta), bam ]
    bam_transcript           = RSEM_CALCULATEEXPRESSION.out.bam_transcript    // channel: [ val(meta), bam ]
    version                  = RSEM_CALCULATEEXPRESSION.out.version           //    path: *.version.txt

    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.csv
    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.csv
    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.csv
    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.csv
}
