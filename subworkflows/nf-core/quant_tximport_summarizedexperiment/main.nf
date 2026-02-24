//
// Quantification post-processing with tximport and SummarizedExperiment
//

include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_UNIFIED       } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT_UNIFIED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

workflow QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT {
    take:
    samplesheet           // channel: [ val(meta), /path/to/samplesheet ]
    quant_results         // channel: [ val(meta), /path/to/results ] - per-sample quant files
    gtf                   // channel: /path/to/genome.gtf
    gtf_id_attribute      //     val: GTF gene ID attribute
    gtf_extra_attribute   //     val: GTF alternative gene attribute (e.g. gene_name)
    quant_type            //     val: 'salmon', 'kallisto', or 'rsem'

    main:
    ch_versions = channel.empty()

    //
    // Create tx2gene mapping from GTF + quantification files
    //
    CUSTOM_TX2GENE (
        gtf.map { gtf_file -> [ [:], gtf_file ] },
        quant_results.collect{ meta_results -> meta_results[1] }.map { results -> [ [:], results ] },
        quant_type,
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    //
    // Import and summarize quantifications with tximport
    //
    TXIMETA_TXIMPORT (
        quant_results.collect{ meta_results -> meta_results[1] }.map { results -> [ ['id': 'all_samples'], results ] },
        CUSTOM_TX2GENE.out.tx2gene,
        quant_type
    )
    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    //
    // Build gene-level SummarizedExperiment
    //
    ch_gene_unified = TXIMETA_TXIMPORT.out.counts_gene
        .join(TXIMETA_TXIMPORT.out.counts_gene_length_scaled, failOnMismatch: true, failOnDuplicate: true)
        .join(TXIMETA_TXIMPORT.out.counts_gene_scaled, failOnMismatch: true, failOnDuplicate: true)
        .join(TXIMETA_TXIMPORT.out.lengths_gene, failOnMismatch: true, failOnDuplicate: true)
        .join(TXIMETA_TXIMPORT.out.tpm_gene, failOnMismatch: true, failOnDuplicate: true)
        .map { row -> tuple(row[0], row.tail()) }

    SE_GENE_UNIFIED (
        ch_gene_unified,
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )
    ch_versions = ch_versions.mix(SE_GENE_UNIFIED.out.versions)

    //
    // Build transcript-level SummarizedExperiment
    //
    ch_transcript_unified = TXIMETA_TXIMPORT.out.counts_transcript
        .join(TXIMETA_TXIMPORT.out.lengths_transcript, failOnMismatch: true, failOnDuplicate: true)
        .join(TXIMETA_TXIMPORT.out.tpm_transcript, failOnMismatch: true, failOnDuplicate: true)
        .map { row -> tuple(row[0], row.tail()) }

    SE_TRANSCRIPT_UNIFIED (
        ch_transcript_unified,
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )
    ch_versions = ch_versions.mix(SE_TRANSCRIPT_UNIFIED.out.versions)

    emit:
    tx2gene                   = CUSTOM_TX2GENE.out.tx2gene                     // channel: [ val(meta), tx2gene.tsv ]

    tpm_gene                  = TXIMETA_TXIMPORT.out.tpm_gene                  //    path: *gene_tpm.tsv
    counts_gene               = TXIMETA_TXIMPORT.out.counts_gene               //    path: *gene_counts.tsv
    lengths_gene              = TXIMETA_TXIMPORT.out.lengths_gene              //    path: *gene_lengths.tsv
    counts_gene_length_scaled = TXIMETA_TXIMPORT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled        = TXIMETA_TXIMPORT.out.counts_gene_scaled        //    path: *gene_counts_scaled.tsv
    tpm_transcript            = TXIMETA_TXIMPORT.out.tpm_transcript            //    path: *transcript_tpm.tsv
    counts_transcript         = TXIMETA_TXIMPORT.out.counts_transcript         //    path: *transcript_counts.tsv
    lengths_transcript        = TXIMETA_TXIMPORT.out.lengths_transcript        //    path: *transcript_lengths.tsv

    merged_gene_rds           = SE_GENE_UNIFIED.out.rds                                //    path: *.rds
    merged_transcript_rds     = SE_TRANSCRIPT_UNIFIED.out.rds                          //    path: *.rds

    versions                  = ch_versions                                    // channel: [ versions.yml ]
}
