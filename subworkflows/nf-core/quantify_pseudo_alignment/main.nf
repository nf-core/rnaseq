//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'
include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE               } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_LENGTH_SCALED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_SCALED        } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT         } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

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
    ch_versions = Channel.empty()

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
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
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
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
    }

    CUSTOM_TX2GENE (
        gtf.map { [ [:], it ] },
        ch_pseudo_results.collect{ it[1] }.map { [ [:], it ] },
        pseudo_aligner,
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT (
        ch_pseudo_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] },
        CUSTOM_TX2GENE.out.tx2gene,
        pseudo_aligner
    )
    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    SE_GENE (
        TXIMETA_TXIMPORT.out.counts_gene.concat(TXIMETA_TXIMPORT.out.tpm_gene).groupTuple(),
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )
    ch_versions = ch_versions.mix(SE_GENE.out.versions)

    SE_GENE_LENGTH_SCALED (
        TXIMETA_TXIMPORT.out.counts_gene_length_scaled.concat(TXIMETA_TXIMPORT.out.tpm_gene).groupTuple(),
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )

    SE_GENE_SCALED (
        TXIMETA_TXIMPORT.out.counts_gene_scaled.concat(TXIMETA_TXIMPORT.out.tpm_gene).groupTuple(),
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )

    SE_TRANSCRIPT (
        TXIMETA_TXIMPORT.out.counts_transcript.concat(TXIMETA_TXIMPORT.out.tpm_transcript).groupTuple(),
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )

    emit:
    results                       = ch_pseudo_results                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                              // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMETA_TXIMPORT.out.tpm_gene                  //    path: *gene_tpm.tsv
    counts_gene                   = TXIMETA_TXIMPORT.out.counts_gene               //    path: *gene_counts.tsv
    lengths_gene                  = TXIMETA_TXIMPORT.out.lengths_gene              //    path: *gene_lengths.tsv
    counts_gene_length_scaled     = TXIMETA_TXIMPORT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled            = TXIMETA_TXIMPORT.out.counts_gene_scaled        //    path: *gene_counts_scaled.tsv
    tpm_transcript                = TXIMETA_TXIMPORT.out.tpm_transcript            //    path: *gene_tpm.tsv
    counts_transcript             = TXIMETA_TXIMPORT.out.counts_transcript         //    path: *transcript_counts.tsv
    lengths_transcript            = TXIMETA_TXIMPORT.out.lengths_transcript        //    path: *transcript_lengths.tsv

    merged_gene_rds               = SE_GENE.out.rds                                //    path: *.rds
    merged_gene_rds_length_scaled = SE_GENE_LENGTH_SCALED.out.rds                  //    path: *.rds
    merged_gene_rds_scaled        = SE_GENE_SCALED.out.rds                         //    path: *.rds
    merged_transcript_rds         = SE_TRANSCRIPT.out.rds                          //    path: *.rds

    versions                      = ch_versions                                    // channel: [ versions.yml ]
}
