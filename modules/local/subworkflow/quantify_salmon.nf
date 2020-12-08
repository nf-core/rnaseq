/*
 * Pseudo-alignment and quantification with Salmon
 */

params.salmon_quant_options = [:]
params.merge_counts_options = [:]

include { SALMON_QUANT        } from '../../nf-core/software/salmon/quant/main' addParams( options: params.salmon_quant_options )
include { SALMON_TX2GENE      } from '../process/salmon_tx2gene'                addParams( options: params.genome_options       )
include { SALMON_TXIMPORT     } from '../process/salmon_tximport'               addParams( options: [publish_by_id : true]      )
include { SALMON_MERGE_COUNTS } from '../process/salmon_merge_counts'           addParams( options: params.merge_counts_options )
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_LENGTH_SCALED
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_SCALED
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_TRANSCRIPT } from '../process/salmon_summarizedexperiment' addParams( options: params.merge_counts_options )

workflow QUANTIFY_SALMON {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    index            // channel: /path/to/salmon/index/
    transcript_fasta // channel: /path/to/transcript.fasta
    gtf              // channel: /path/to/genome.gtf
    alignment_mode   //    bool: Run Salmon in alignment mode

    main:    
    /*
     * Quantify and merge counts across samples
     */
    SALMON_QUANT        ( reads, index, gtf, transcript_fasta, alignment_mode)
    SALMON_TX2GENE      ( SALMON_QUANT.out.results.collect{it[1]}, gtf )
    SALMON_TXIMPORT     ( SALMON_QUANT.out.results, SALMON_TX2GENE.out.collect() )
    SALMON_MERGE_COUNTS (
        SALMON_TXIMPORT.out.counts_gene.collect{it[1]},                // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene.collect{it[1]},
        SALMON_TXIMPORT.out.counts_gene_length_scaled.collect{it[1]},  // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene_length_scaled.collect{it[1]},
        SALMON_TXIMPORT.out.counts_gene_scaled.collect{it[1]},         // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene_scaled.collect{it[1]},
        SALMON_TXIMPORT.out.counts_transcript.collect{it[1]},   
        SALMON_TXIMPORT.out.tpm_transcript.collect{it[1]},
    )

    SALMON_SE_GENE ( 
        SALMON_MERGE_COUNTS.out.counts_gene,
        SALMON_MERGE_COUNTS.out.tpm_gene,
        SALMON_TX2GENE.out.collect(),
    )

    SALMON_SE_GENE_LENGTH_SCALED ( 
        SALMON_MERGE_COUNTS.out.counts_gene_length_scaled,
        SALMON_MERGE_COUNTS.out.tpm_gene_length_scaled,
        SALMON_TX2GENE.out.collect(),
    )

    SALMON_SE_GENE_SCALED ( 
        SALMON_MERGE_COUNTS.out.counts_gene_scaled,
        SALMON_MERGE_COUNTS.out.tpm_gene_scaled,
        SALMON_TX2GENE.out.collect(),
    )

    SALMON_SE_TRANSCRIPT ( 
        SALMON_MERGE_COUNTS.out.counts_transcript,
        SALMON_MERGE_COUNTS.out.tpm_transcript,
        SALMON_TX2GENE.out.collect(),
    )

    emit:
    results                          = SALMON_QUANT.out.results                             // channel: [ val(meta), results_dir ]
    salmon_version                   = SALMON_QUANT.out.version                             //    path: *.version.txt
    
    tpm_gene                         = SALMON_TXIMPORT.out.tpm_gene                         // channel: [ val(meta), counts ]
    counts_gene                      = SALMON_TXIMPORT.out.counts_gene                      // channel: [ val(meta), counts ]
    tpm_gene_length_scaled           = SALMON_TXIMPORT.out.tpm_gene_length_scaled           // channel: [ val(meta), counts ]
    counts_gene_length_scaled        = SALMON_TXIMPORT.out.counts_gene_length_scaled        // channel: [ val(meta), counts ]
    tpm_gene_scaled                  = SALMON_TXIMPORT.out.tpm_gene_scaled                  // channel: [ val(meta), counts ]
    counts_gene_scaled               = SALMON_TXIMPORT.out.counts_gene_scaled               // channel: [ val(meta), counts ]
    tpm_transcript                   = SALMON_TXIMPORT.out.tpm_transcript                   // channel: [ val(meta), counts ]
    counts_transcript                = SALMON_TXIMPORT.out.counts_transcript                // channel: [ val(meta), counts ]
    tximeta_version                  = SALMON_TXIMPORT.out.version                          //    path: *.version.txt
    
    merged_counts_gene               = SALMON_MERGE_COUNTS.out.counts_gene                  //    path: *.gene_counts.tsv
    merged_tpm_gene                  = SALMON_MERGE_COUNTS.out.tpm_gene                     //    path: *.gene_tpm.tsv
    merged_counts_gene_length_scaled = SALMON_MERGE_COUNTS.out.counts_gene_length_scaled    //    path: *.gene_counts.tsv
    merged_tpm_gene_length_scaled    = SALMON_MERGE_COUNTS.out.tpm_gene_length_scaled       //    path: *.gene_tpm.tsv
    merged_counts_gene_scaled        = SALMON_MERGE_COUNTS.out.counts_gene_scaled           //    path: *.gene_counts.tsv
    merged_tpm_gene_scaled           = SALMON_MERGE_COUNTS.out.tpm_gene_scaled              //    path: *.gene_tpm.tsv
    merged_gene_rds                  = SALMON_SE_GENE.out.rds                               //    path: *.rds
    merged_gene_rds_length_scaled    = SALMON_SE_GENE_LENGTH_SCALED.out.rds                 //    path: *.rds
    merged_gene_rds_scaled           = SALMON_SE_GENE_SCALED.out.rds                        //    path: *.rds
    summarizedexperiment_version     = SALMON_SE_GENE.out.version                           //    path: *.version.txt
        
    merged_counts_transcript         = SALMON_MERGE_COUNTS.out.counts_transcript            //    path: *.transcript_counts.tsv
    merged_tpm_transcript            = SALMON_MERGE_COUNTS.out.tpm_transcript               //    path: *.transcript_tpm.tsv
    merged_transcript_rds            = SALMON_SE_TRANSCRIPT.out.rds                         //    path: *.rds
}
