/*
 * Pseudo-alignment and quantification with Salmon
 */

params.index_options        = [:]
params.genome_options       = [:]
params.salmon_index_options = [:]
params.salmon_quant_options = [:]
params.merge_counts_options = [:]

include { GUNZIP              } from '../process/gunzip'                        addParams( options: params.genome_options       )
include { UNTAR               } from '../process/untar'                         addParams( options: params.index_options        )
include { GTF_GENE_FILTER     } from '../process/gtf_gene_filter'               addParams( options: params.genome_options       )
include { GFFREAD             } from '../process/gffread'                       addParams( options: params.genome_options       )
include { SALMON_INDEX        } from '../../nf-core/software/salmon/index/main' addParams( options: params.salmon_index_options )
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
    index            //    file: /path/to/salmon/index/
    transcript_fasta //    file: /path/to/transcript.fasta
    genome_fasta     //    file: /path/to/genome.fasta
    gtf              //    file: /path/to/genome.gtf

    main:
    /*
     * Uncompress Salmon index or generate from scratch if required
     */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        if (transcript_fasta) {
            if (transcript_fasta.endsWith('.gz')) {
                ch_transcript_fasta = GUNZIP ( transcript_fasta ).gunzip
            } else {
                ch_transcript_fasta = file(transcript_fasta)
            }
        } else {
            ch_transcript_fasta = GFFREAD ( genome_fasta, GTF_GENE_FILTER ( genome_fasta, gtf ) ).fasta
        }
        ch_index = SALMON_INDEX ( ch_transcript_fasta ).index
    }

    /*
     * Quantify and merge counts across samples
     */
    SALMON_QUANT        ( reads, ch_index, gtf )
    SALMON_TX2GENE      ( SALMON_QUANT.out.results.collect{it[1]}, gtf )
    SALMON_TXIMPORT     ( SALMON_QUANT.out.results, SALMON_TX2GENE.out.collect() )
    SALMON_MERGE_COUNTS (
        SALMON_TXIMPORT.out.counts_gene.collect{it[1]},                       // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene.collect{it[1]},
        SALMON_TXIMPORT.out.counts_gene_length_scaled.collect{it[1]},         // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene_length_scaled.collect{it[1]},
        SALMON_TXIMPORT.out.counts_gene_scaled.collect{it[1]},                // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        SALMON_TXIMPORT.out.tpm_gene_scaled.collect{it[1]},
        SALMON_TXIMPORT.out.counts_transcript.collect{it[1]},   
        SALMON_TXIMPORT.out.tpm_transcript.collect{it[1]},
    )
    SALMON_SE_GENE ( 
        SALMON_MERGE_COUNTS.out.counts_gene,
        SALMON_MERGE_COUNTS.out.tpm_gene,
        SALMON_TX2GENE.out.collect(),
    )
    SALMON_SE_GENE_LS ( 
        SALMON_MERGE_COUNTS.out.counts_gene_length_scaled,
        SALMON_MERGE_COUNTS.out.tpm_gene_length_scaled,
        SALMON_TX2GENE.out.collect(),
    )
    SALMON_SE_GENE_S ( 
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
    results                      = SALMON_QUANT.out.results                             // channel: [ val(meta), results_dir ]
    salmon_version               = SALMON_QUANT.out.version                             //    path: *.version.txt
    
    tpm_gene                     = SALMON_TXIMPORT.out.tpm_gene                         // channel: [ val(meta), counts ]
    counts_gene                  = SALMON_TXIMPORT.out.counts_gene                      // channel: [ val(meta), counts ]
    tpm_gene_ls                  = SALMON_TXIMPORT.out.tpm_gene_length_scaled           // channel: [ val(meta), counts ]
    counts_gene_ls               = SALMON_TXIMPORT.out.counts_gene_length_scaled        // channel: [ val(meta), counts ]
    tpm_gene_s                   = SALMON_TXIMPORT.out.tpm_gene_scaled                  // channel: [ val(meta), counts ]
    counts_gene_s                = SALMON_TXIMPORT.out.counts_gene_scaled               // channel: [ val(meta), counts ]
    tpm_transcript               = SALMON_TXIMPORT.out.tpm_transcript                   // channel: [ val(meta), counts ]
    counts_transcript            = SALMON_TXIMPORT.out.counts_transcript                // channel: [ val(meta), counts ]
    tximeta_version              = SALMON_TXIMPORT.out.version                          //    path: *.version.txt
    
    merged_counts_gene           = SALMON_MERGE_COUNTS.out.counts_gene                  //    path: *.gene_counts.tsv
    merged_tpm_gene              = SALMON_MERGE_COUNTS.out.tpm_gene                     //    path: *.gene_tpm.tsv
    merged_counts_gene_ls        = SALMON_MERGE_COUNTS.out.counts_gene_length_scaled    //    path: *.gene_counts.tsv
    merged_tpm_gene_s            = SALMON_MERGE_COUNTS.out.tpm_gene_length_scaled       //    path: *.gene_tpm.tsv
    merged_counts_gene_ls        = SALMON_MERGE_COUNTS.out.counts_gene_scaled           //    path: *.gene_counts.tsv
    merged_tpm_gene_s            = SALMON_MERGE_COUNTS.out.tpm_gene_scaled              //    path: *.gene_tpm.tsv
    merged_gene_rds              = SALMON_SE_GENE.out.rds                               //    path: *.rds
    merged_gene_rds_ls           = SALMON_SE_GENE_LENGTH_SCALED.out.rds                 //    path: *.rds
    merged_gene_rds_s            = SALMON_SE_GENE_SCALED.out.rds                        //    path: *.rds
    summarizedexperiment_version = SALMON_SE_GENE.out.version                           //    path: *.version.txt
        
    merged_counts_transcript     = SALMON_MERGE_COUNTS.out.counts_transcript            //    path: *.transcript_counts.tsv
    merged_tpm_transcript        = SALMON_MERGE_COUNTS.out.tpm_transcript               //    path: *.transcript_tpm.tsv
    merged_transcript_rds        = SALMON_SE_TRANSCRIPT.out.rds                         //    path: *.rds
}
