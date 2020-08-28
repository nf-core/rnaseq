/*
 * Pseudo-alignment and quantification with Salmon
 */

include { GUNZIP            } from '../process/gunzip'
include { UNTAR             } from '../process/untar'
include { TRANSCRIPTS2FASTA } from '../process/transcripts2fasta'
include { SALMON_INDEX      } from '../process/salmon_index'
include { SALMON_QUANT      } from '../process/salmon_quant'
include { SALMON_TX2GENE    } from '../process/salmon_tx2gene'
include { SALMON_TXIMPORT   } from '../process/salmon_tximport'
include { SALMON_MERGE      } from '../process/salmon_merge'

workflow QUANTIFY_SALMON {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index             // file: /path/to/salmon/index/
    fasta             // file: /path/to/transcripts.fasta
    gtf               // file: /path/to/genome.gtf
    prep_options     //  map: module options for additional arguments and publishing files
    index_options
    quant_options

    main:
    /*
     * Uncompress Salmon index or generate from scratch if required
     */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, prep_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        if (fasta) {
            if (fasta.endsWith('.gz')) {
                ch_fasta = GUNZIP ( fasta, prep_options ).gunzip
            } else {
                ch_fasta = file(fasta)
            }
        } else {
            ch_fasta = TRANSCRIPTS2FASTA ( fasta, gtf, prep_options ).fasta
        }
        ch_index = SALMON_INDEX ( ch_fasta, index_options )
    }

    /*
     * Quantify and merge counts across samples
     */
    SALMON_QUANT    ( reads, ch_index, gtf, quant_options )
    SALMON_TX2GENE  ( SALMON_QUANT.out.results.collect{it[1]}, ch_gtf, prep_options )
    SALMON_TXIMPORT ( SALMON_QUANT.out.results, SALMON_TX2GENE.out.collect(), [publish_by_id : true] )
    SALMON_MERGE    (
        SALMON_TXIMPORT.out.gene_tpm.collect{it[1]},
        SALMON_TXIMPORT.out.gene_counts.collect{it[1]},
        SALMON_TXIMPORT.out.transcript_tpm.collect{it[1]},
        SALMON_TXIMPORT.out.transcript_counts.collect{it[1]},
        SALMON_TX2GENE.out.collect(),
        [:]
    )

    emit:
    results        = SALMON_QUANT.out.results // channel: [ val(meta), results_dir ]
    salmon_version = SALMON_QUANT.out.version // path: *.version.txt
}
