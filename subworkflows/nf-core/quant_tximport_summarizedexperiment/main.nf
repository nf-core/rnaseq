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
    skip_merge            //    bool: skip cross-sample merging, run tximport per-sample

    main:

    //
    // Create tx2gene mapping from GTF + a single sample's quantification files.
    // tx2gene only uses quant files to discover which GTF attribute corresponds
    // to transcript IDs, so a single sample suffices. This assumes all samples
    // were quantified against the same transcriptome, avoiding the need to run
    // tx2gene independently per sample. If a use case arises requiring mixed
    // transcriptomes, this will need to be revisited.
    //
    ch_tx2gene_quants = quant_results
        .first()
        .map { meta, results -> [ [:], results ] }

    CUSTOM_TX2GENE (
        gtf.map { gtf_file -> [ [:], gtf_file ] },
        ch_tx2gene_quants,
        quant_type,
        gtf_id_attribute,
        gtf_extra_attribute
    )

    //
    // Import and summarize quantifications with tximport
    // In per-sample mode, run once per sample instead of collecting all
    //
    ch_tximport_input = skip_merge
        ? quant_results
        : quant_results.collect{ meta_results -> meta_results[1] }.map { results -> [ ['id': 'all_samples'], results ] }

    TXIMETA_TXIMPORT (
        ch_tximport_input,
        CUSTOM_TX2GENE.out.tx2gene,
        quant_type
    )

    //
    // Build SummarizedExperiment objects (only when merging)
    //
    ch_merged_gene_rds       = channel.empty()
    ch_merged_transcript_rds = channel.empty()

    if (!skip_merge) {
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

        ch_merged_gene_rds       = SE_GENE_UNIFIED.out.rds
        ch_merged_transcript_rds = SE_TRANSCRIPT_UNIFIED.out.rds
    }

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

    merged_gene_rds           = ch_merged_gene_rds                             //    path: *.rds
    merged_transcript_rds     = ch_merged_transcript_rds                       //    path: *.rds
}
