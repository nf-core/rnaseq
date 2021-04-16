/*
 * Pseudo-alignment and quantification with Salmon
 */

params.genome_options       = [:]
params.salmon_quant_options = [:]
params.tximport_options     = [:]
params.merge_counts_options = [:]

include { SALMON_QUANT        } from '../../modules/nf-core/software/salmon/quant/main' addParams( options: params.salmon_quant_options )
include { SALMON_TX2GENE      } from '../../modules/local/salmon_tx2gene'               addParams( options: params.genome_options       )
include { SALMON_TXIMPORT     } from '../../modules/local/salmon_tximport'              addParams( options: params.tximport_options     )
include { SALMON_MERGE_COUNTS } from '../../modules/local/salmon_merge_counts'          addParams( options: params.merge_counts_options )
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_LENGTH_SCALED
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_SCALED
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_TRANSCRIPT } from '../../modules/local/salmon_summarizedexperiment' addParams( options: params.merge_counts_options )

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
    SALMON_QUANT        ( reads, index, gtf, transcript_fasta, alignment_mode )
    SALMON_TX2GENE      ( SALMON_QUANT.out.results.collect{it[1]}, gtf )
    SALMON_TXIMPORT     ( SALMON_QUANT.out.results.collect{it[1]}, SALMON_TX2GENE.out.collect() )


    SALMON_SE_GENE ( 
        SALMON_TXIMPORT.out.counts_gene,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.collect()
    )

    SALMON_SE_GENE_LENGTH_SCALED ( 
        SALMON_TXIMPORT.out.counts_gene_length_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.collect()
    )

    SALMON_SE_GENE_SCALED ( 
        SALMON_TXIMPORT.out.counts_gene_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.collect()
    )

    SALMON_SE_TRANSCRIPT ( 
        SALMON_TXIMPORT.out.counts_transcript,
        SALMON_TXIMPORT.out.tpm_transcript,
        SALMON_TX2GENE.out.collect()
    )

    emit:
    results                          = SALMON_QUANT.out.results                      // channel: [ val(meta), results_dir ]
    salmon_version                   = SALMON_QUANT.out.version                      //    path: *.version.txt
    
    tpm_gene                         = SALMON_TXIMPORT.out.tpm_gene                  // channel: [ val(meta), counts ]
    counts_gene                      = SALMON_TXIMPORT.out.counts_gene               // channel: [ val(meta), counts ]
    counts_gene_length_scaled        = SALMON_TXIMPORT.out.counts_gene_length_scaled // channel: [ val(meta), counts ]
    counts_gene_scaled               = SALMON_TXIMPORT.out.counts_gene_scaled        // channel: [ val(meta), counts ]
    tpm_transcript                   = SALMON_TXIMPORT.out.tpm_transcript            // channel: [ val(meta), counts ]
    counts_transcript                = SALMON_TXIMPORT.out.counts_transcript         // channel: [ val(meta), counts ]
    tximeta_version                  = SALMON_TXIMPORT.out.version                   //    path: *.version.txt
    
    merged_gene_rds                  = SALMON_SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled    = SALMON_SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled           = SALMON_SE_GENE_SCALED.out.rds                 //    path: *.rds
    summarizedexperiment_version     = SALMON_SE_GENE.out.version                    //    path: *.version.txt
        
    merged_counts_transcript         = SALMON_TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript            = SALMON_TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds            = SALMON_SE_TRANSCRIPT.out.rds                  //    path: *.rds
}
