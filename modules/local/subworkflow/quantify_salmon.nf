/*
 * Pseudo-alignment and quantification with Salmon
 */

include { GUNZIP              } from '../process/gunzip'
include { UNTAR               } from '../process/untar'
include { GTF_GENE_FILTER     } from '../process/gtf_gene_filter'
include { GFFREAD             } from '../process/gffread'
include { SALMON_TX2GENE      } from '../process/salmon_tx2gene'
include { SALMON_TXIMPORT     } from '../process/salmon_tximport'
include { SALMON_MERGE_COUNTS } from '../process/salmon_merge_counts'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE 
          SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_TRANSCRIPT } from '../process/salmon_summarizedexperiment'

include { SALMON_INDEX        } from '../../nf-core/software/salmon/index/main'
include { SALMON_QUANT        } from '../../nf-core/software/salmon/quant/main'

workflow QUANTIFY_SALMON {
    take:
    reads                  // channel: [ val(meta), [ reads ] ]
    index                  //    file: /path/to/salmon/index/
    transcript_fasta       //    file: /path/to/transcript.fasta
    genome_fasta           //    file: /path/to/genome.fasta
    gtf                    //    file: /path/to/genome.gtf
    publish_index_options  //     map: options for publishing index
    publish_genome_options //     map: options for publishing genome files
    salmon_index_options   //     map: options for salmon_index module
    salmon_quant_options   //     map: options for salmon_quant module
    merge_counts_options   //     map: options for merge_counts_salmon module

    main:
    /*
     * Uncompress Salmon index or generate from scratch if required
     */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, publish_index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        if (transcript_fasta) {
            if (transcript_fasta.endsWith('.gz')) {
                ch_transcript_fasta = GUNZIP ( transcript_fasta, publish_genome_options ).gunzip
            } else {
                ch_transcript_fasta = file(transcript_fasta)
            }
        } else {
            ch_transcript_fasta = GFFREAD ( genome_fasta, GTF_GENE_FILTER ( genome_fasta, gtf, publish_genome_options ), publish_genome_options).fasta
        }
        ch_index = SALMON_INDEX ( ch_transcript_fasta, salmon_index_options ).index
    }

    /*
     * Quantify and merge counts across samples
     */
    SALMON_QUANT        ( reads, ch_index, gtf, salmon_quant_options )
    SALMON_TX2GENE      ( SALMON_QUANT.out.results.collect{it[1]}, gtf, publish_genome_options )
    SALMON_TXIMPORT     ( SALMON_QUANT.out.results, SALMON_TX2GENE.out.collect(), [publish_by_id : true] )
    SALMON_MERGE_COUNTS (
        SALMON_TXIMPORT.out.counts_gene.collect{it[1]},
        SALMON_TXIMPORT.out.tpm_gene.collect{it[1]},
        SALMON_TXIMPORT.out.counts_transcript.collect{it[1]},
        SALMON_TXIMPORT.out.tpm_transcript.collect{it[1]},
        merge_counts_options
    )
    SALMON_SE_GENE ( 
        SALMON_MERGE_COUNTS.out.counts_gene,
        SALMON_MERGE_COUNTS.out.tpm_gene,
        SALMON_TX2GENE.out.collect(),
        merge_counts_options
    )
    SALMON_SE_TRANSCRIPT ( 
        SALMON_MERGE_COUNTS.out.counts_transcript,
        SALMON_MERGE_COUNTS.out.tpm_transcript,
        SALMON_TX2GENE.out.collect(),
        merge_counts_options
    )

    emit:
    results                      = SALMON_QUANT.out.results                  // channel: [ val(meta), results_dir ]
    salmon_version               = SALMON_QUANT.out.version                  //    path: *.version.txt

    tpm_gene                     = SALMON_TXIMPORT.out.tpm_gene              // channel: [ val(meta), counts ]
    counts_gene                  = SALMON_TXIMPORT.out.counts_gene           // channel: [ val(meta), counts ]
    tpm_transcript               = SALMON_TXIMPORT.out.tpm_transcript        // channel: [ val(meta), counts ]
    counts_transcript            = SALMON_TXIMPORT.out.counts_transcript     // channel: [ val(meta), counts ]
    tximeta_version              = SALMON_TXIMPORT.out.version               //    path: *.version.txt

    merged_counts_gene           = SALMON_MERGE_COUNTS.out.counts_gene       //    path: *.gene_counts.tsv
    merged_tpm_gene              = SALMON_MERGE_COUNTS.out.tpm_gene          //    path: *.gene_tpm.tsv
    merged_gene_rds              = SALMON_SE_GENE.out.rds                    //    path: *.rds
    summarizedexperiment_version = SALMON_SE_GENE.out.version                //    path: *.version.txt
    
    merged_counts_transcript     = SALMON_MERGE_COUNTS.out.counts_transcript //    path: *.transcript_counts.tsv
    merged_tpm_transcript        = SALMON_MERGE_COUNTS.out.tpm_transcript    //    path: *.transcript_tpm.tsv
    merged_transcript_rds        = SALMON_SE_TRANSCRIPT.out.rds              //    path: *.rds
}
