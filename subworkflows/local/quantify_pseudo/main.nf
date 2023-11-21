//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT   } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT } from '../../../modules/nf-core/kallisto/quant'
include { TX2GENE        } from '../../../modules/local/tx2gene'
include { TXIMPORT       } from '../../../modules/local/tximport'

include { SUMMARIZEDEXPERIMENT as SE_GENE               } from '../../../modules/local/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT as SE_GENE_LENGTH_SCALED } from '../../../modules/local/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT as SE_GENE_SCALED        } from '../../../modules/local/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT         } from '../../../modules/local/summarizedexperiment'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: /path/to//index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    pseudo_aligner            //     val: kallisto or salmon 
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode 

    main:

    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT ( reads, index, gtf, transcript_fasta, alignment_mode, lib_type )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
    } else {
        KALLISTO_QUANT ( reads, index, gtf, [], kallisto_quant_fraglen, kallisto_quant_fraglen_sd)
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
    }

    TX2GENE ( ch_pseudo_results.collect{it[1]}, pseudo_aligner, gtf )
    ch_versions = ch_versions.mix(TX2GENE.out.versions)

    TXIMPORT ( ch_pseudo_results.collect{it[1]}, TX2GENE.out.tsv.collect(), pseudo_aligner )
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    SE_GENE (
        TXIMPORT.out.counts_gene,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )
    ch_versions = ch_versions.mix(SE_GENE.out.versions)

    SE_GENE_LENGTH_SCALED (
        TXIMPORT.out.counts_gene_length_scaled,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )

    SE_GENE_SCALED (
        TXIMPORT.out.counts_gene_scaled,
        TXIMPORT.out.tpm_gene,
        TX2GENE.out.tsv.collect()
    )

    SE_TRANSCRIPT (
        TXIMPORT.out.counts_transcript,
        TXIMPORT.out.tpm_transcript,
        TX2GENE.out.tsv.collect()
    )

    emit:
    results                       = ch_pseudo_results                      // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                      // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMPORT.out.tpm_gene                  //    path *gene_tpm.tsv
    counts_gene                   = TXIMPORT.out.counts_gene               //    path *gene_counts.tsv
    lengths_gene                  = TXIMPORT.out.lengths_gene              //    path *gene_lengths.tsv
    counts_gene_length_scaled     = TXIMPORT.out.counts_gene_length_scaled //    path *gene_counts_length_scaled.tsv
    counts_gene_scaled            = TXIMPORT.out.counts_gene_scaled        //    path *gene_counts_scaled.tsv
    tpm_transcript                = TXIMPORT.out.tpm_transcript            //    path *gene_tpm.tsv
    counts_transcript             = TXIMPORT.out.counts_transcript         //    path *transcript_counts.tsv
    lengths_transcript            = TXIMPORT.out.lengths_transcript        //    path *transcript_lengths.tsv

    merged_gene_rds               = SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled = SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled        = SE_GENE_SCALED.out.rds                 //    path: *.rds

    merged_counts_transcript      = TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds         = SE_TRANSCRIPT.out.rds                  //    path: *.rds

    versions                      = ch_versions                            // channel: [ versions.yml ]
}
