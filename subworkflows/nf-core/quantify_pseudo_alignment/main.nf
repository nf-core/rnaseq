//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: /path/to//index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner            //     val: kallisto or salmon
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode

    main:
    ch_versions = channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT (
            reads,
            index,
            gtf,
            transcript_fasta,
            alignment_mode,
            lib_type
        )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
    } else {
        KALLISTO_QUANT (
            reads,
            index,
            gtf,
            [],
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT (
        samplesheet,
        ch_pseudo_results,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner
    )
    ch_versions = ch_versions.mix(QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.versions)

    emit:
    results                       = ch_pseudo_results                                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                                              // channel: [ val(meta), files_for_multiqc ]
    tx2gene                       = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tx2gene                // channel: [ val(meta), tx2gene.tsv ]

    tpm_gene                      = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_gene               //    path: *gene_tpm.tsv
    counts_gene                   = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene            //    path: *gene_counts.tsv
    lengths_gene                  = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_gene           //    path: *gene_lengths.tsv
    counts_gene_length_scaled     = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled            = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_scaled     //    path: *gene_counts_scaled.tsv
    tpm_transcript                = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_transcript         //    path: *transcript_tpm.tsv
    counts_transcript             = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_transcript      //    path: *transcript_counts.tsv
    lengths_transcript            = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_transcript     //    path: *transcript_lengths.tsv

    merged_gene_rds_unified       = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_gene_rds       //    path: *.rds
    merged_transcript_rds_unified = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_transcript_rds //    path: *.rds

    versions                      = ch_versions                                                    // channel: [ versions.yml ]
}
