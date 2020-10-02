/*
 * Pseudo-alignment and quantification with Salmon
 */

include { GUNZIP              } from '../process/gunzip'
include { UNTAR               } from '../process/untar'
include { TRANSCRIPTS2FASTA   } from '../process/transcripts2fasta'
include { SALMON_TX2GENE      } from '../process/salmon_tx2gene'
include { SALMON_TXIMPORT     } from '../process/salmon_tximport'
include { SALMON_MERGE_COUNTS } from '../process/salmon_merge_counts'
include { SALMON_INDEX        } from '../../nf-core/software/salmon/index/main'
include { SALMON_QUANT        } from '../../nf-core/software/salmon/quant/main'

workflow QUANTIFY_SALMON {
    take:
    reads                  // channel: [ val(meta), [ reads ] ]
    index                  //    file: /path/to/salmon/index/
    fasta                  //    file: /path/to/transcripts.fasta
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
        if (fasta) {
            if (fasta.endsWith('.gz')) {
                ch_fasta = GUNZIP ( fasta, publish_genome_options ).gunzip
            } else {
                ch_fasta = file(fasta)
            }
        } else {
            ch_fasta = TRANSCRIPTS2FASTA ( fasta, gtf, publish_genome_options ).fasta
        }
        ch_index = SALMON_INDEX ( ch_fasta, salmon_index_options ).index
    }

    /*
     * Quantify and merge counts across samples
     */
    SALMON_QUANT        ( reads, ch_index, gtf, salmon_quant_options )
    SALMON_TX2GENE      ( SALMON_QUANT.out.results.collect{it[1]}, gtf, publish_genome_options )
    SALMON_TXIMPORT     ( SALMON_QUANT.out.results, SALMON_TX2GENE.out.collect(), [publish_by_id : true] )
    SALMON_MERGE_COUNTS (
        SALMON_TXIMPORT.out.tpm_gene.collect{it[1]},
        SALMON_TXIMPORT.out.counts_gene.collect{it[1]},
        SALMON_TXIMPORT.out.tpm_transcript.collect{it[1]},
        SALMON_TXIMPORT.out.counts_transcript.collect{it[1]},
        SALMON_TX2GENE.out.collect(),
        merge_counts_options
    )

    emit:
    results                  = SALMON_QUANT.out.results                  // channel: [ val(meta), results_dir ]
    version                  = SALMON_QUANT.out.version                  //    path: *.version.txt

    tpm_gene                 = SALMON_TXIMPORT.out.tpm_gene              // channel: [ val(meta), counts ]
    counts_gene              = SALMON_TXIMPORT.out.counts_gene           // channel: [ val(meta), counts ]
    tpm_transcript           = SALMON_TXIMPORT.out.tpm_transcript        // channel: [ val(meta), counts ]
    counts_transcript        = SALMON_TXIMPORT.out.counts_transcript     // channel: [ val(meta), counts ]

    merged_tpm_gene          = SALMON_MERGE_COUNTS.out.tpm_gene          //    path: *.gene_tpm.csv
    merged_counts_gene       = SALMON_MERGE_COUNTS.out.counts_gene       //    path: *.gene_counts.csv
    merged_tpm_transcript    = SALMON_MERGE_COUNTS.out.tpm_transcript    //    path: *.transcript_tpm.csv
    merged_counts_transcript = SALMON_MERGE_COUNTS.out.counts_transcript //    path: *.transcript_counts.csv
    rds                      = SALMON_MERGE_COUNTS.out.rds               //    path: *.rds
}
