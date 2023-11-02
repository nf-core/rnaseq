//
// Pseudo-alignment and quantification with Kallisto
//

include { KALLISTO_QUANT    } from '../../../modules/nf-core/kallisto/quant'
include { TX2GENE  } from '../../../modules/local/tx2gene'
include { TXIMPORT } from '../../../modules/local/tximport'

include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE               } from '../../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_LENGTH_SCALED } from '../../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_SCALED        } from '../../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_TRANSCRIPT         } from '../../../modules/local/salmon_summarizedexperiment'

workflow QUANTIFY_KALLISTO {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    index            // channel: [ val(meta2), /path/to/kallisto/index/ ]
    transcript_fasta // channel: /path/to/transcript.fasta
    gtf              // channel: /path/to/genome.gtf
    alignment_mode   //    bool: Run Salmon in alignment mode
    lib_type         //     val: String to override salmon library type

    main:

    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    reads.view()
    KALLISTO_QUANT ( reads, index, gtf, [])
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())

    TX2GENE ( KALLISTO_QUANT.out.results.collect{it[1]}, 'kallisto', gtf )
    ch_versions = ch_versions.mix(TX2GENE.out.versions)

    TXIMPORT ( KALLISTO_QUANT.out.results.collect{it[1]}, TX2GENE.out.tsv.collect(), 'kallisto' )
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    SALMON_SE_GENE (
        TXIMPORT.out.counts_gene,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )
    ch_versions = ch_versions.mix(SALMON_SE_GENE.out.versions)

    SALMON_SE_GENE_LENGTH_SCALED (
        TXIMPORT.out.counts_gene_length_scaled,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )

    SALMON_SE_GENE_SCALED (
        TXIMPORT.out.counts_gene_scaled,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )

    SALMON_SE_TRANSCRIPT (
        TXIMPORT.out.counts_transcript,
        TXIMPORT.out.tpm_transcript,
        TX2GENE.out.tsv.collect()
    )

    emit:
    results                       = KALLISTO_QUANT.out.results                    // channel: [ val(meta), abundances ]

    tpm_gene                      = TXIMPORT.out.tpm_gene                  // channel: [ val(meta), counts ]
    counts_gene                   = TXIMPORT.out.counts_gene               // channel: [ val(meta), counts ]
    counts_gene_length_scaled     = TXIMPORT.out.counts_gene_length_scaled // channel: [ val(meta), counts ]
    counts_gene_scaled            = TXIMPORT.out.counts_gene_scaled        // channel: [ val(meta), counts ]
    tpm_transcript                = TXIMPORT.out.tpm_transcript            // channel: [ val(meta), counts ]
    counts_transcript             = TXIMPORT.out.counts_transcript         // channel: [ val(meta), counts ]

    merged_gene_rds               = SALMON_SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled = SALMON_SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled        = SALMON_SE_GENE_SCALED.out.rds                 //    path: *.rds

    merged_counts_transcript      = TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds         = SALMON_SE_TRANSCRIPT.out.rds                  //    path: *.rds

    versions                      = ch_versions                                   // channel: [ versions.yml ]
}
