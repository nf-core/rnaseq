//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'

workflow QUANTIFY_RSEM {
    take:
    samplesheet           // channel: [ val(meta), /path/to/samplesheet ]
    reads                 // channel: [ val(meta), [ reads ] ] - FASTQ or BAM files
    index                 // channel: /path/to/rsem/index/
    gtf                   // channel: /path/to/genome.gtf
    gtf_id_attribute      //     val: GTF gene ID attribute
    gtf_extra_attribute   //     val: GTF alternative gene attribute (e.g. gene_name)
    use_sentieon_star     // boolean: use Sentieon-accelerated STAR (FASTQ mode only)

    main:

    ch_versions = channel.empty()

    //
    // Quantify reads with RSEM
    //
    ch_rsem_out = null
    if (use_sentieon_star) {
        SENTIEON_RSEMCALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = SENTIEON_RSEMCALCULATEEXPRESSION
    } else {
        RSEM_CALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = RSEM_CALCULATEEXPRESSION
    }

    ch_counts_gene       = ch_rsem_out.out.counts_gene
    ch_counts_transcript = ch_rsem_out.out.counts_transcript
    ch_stat              = ch_rsem_out.out.stat
    ch_logs              = ch_rsem_out.out.logs

    //
    // Merge counts across samples
    //
    CUSTOM_RSEMMERGECOUNTS (
        ch_counts_gene.collect{ it[1] }.map { results -> [ ['id': 'all_samples'], results ] },
        ch_counts_transcript.collect{ it[1] }
    )

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT (
        samplesheet,
        ch_counts_transcript,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        'rsem'
    )
    ch_versions = ch_versions.mix(QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.versions)

    emit:
    // Per-sample outputs
    raw_counts_gene          = ch_counts_gene                                                        // channel: [ val(meta), counts ]
    raw_counts_transcript    = ch_counts_transcript                                                  // channel: [ val(meta), counts ]
    stat                     = ch_stat                                                               // channel: [ val(meta), stat ]
    logs                     = ch_logs                                                               // channel: [ val(meta), logs ]

    // RSEM merge outputs
    merged_counts_gene       = CUSTOM_RSEMMERGECOUNTS.out.counts_gene                                // channel: [ val(meta), counts ]
    merged_tpm_gene          = CUSTOM_RSEMMERGECOUNTS.out.tpm_gene                                   // channel: [ val(meta), tpm ]
    merged_counts_transcript = CUSTOM_RSEMMERGECOUNTS.out.counts_transcript                          // channel: [ val(meta), counts ]
    merged_tpm_transcript    = CUSTOM_RSEMMERGECOUNTS.out.tpm_transcript                             // channel: [ val(meta), tpm ]
    merged_genes_long        = CUSTOM_RSEMMERGECOUNTS.out.genes_long                                 // channel: [ val(meta), genes_long ]
    merged_isoforms_long     = CUSTOM_RSEMMERGECOUNTS.out.isoforms_long                              // channel: [ val(meta), isoforms_long ]

    // tximport outputs
    tpm_gene                  = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_gene                     //    path: *gene_tpm.tsv
    counts_gene               = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene                  //    path: *gene_counts.tsv
    counts_gene_length_scaled = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_length_scaled    //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled        = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_scaled           //    path: *gene_counts_scaled.tsv
    lengths_gene              = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_gene                 //    path: *gene_lengths.tsv
    tpm_transcript            = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_transcript               //    path: *transcript_tpm.tsv
    counts_transcript         = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_transcript            //    path: *transcript_counts.tsv
    lengths_transcript        = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_transcript           //    path: *transcript_lengths.tsv

    // SummarizedExperiment objects
    merged_gene_rds_unified       = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_gene_rds          //    path: *.rds
    merged_transcript_rds_unified = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_transcript_rds    //    path: *.rds

    // tx2gene
    tx2gene                   = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tx2gene                      //    path: *tx2gene.tsv

    versions                  = ch_versions                                                          // channel: [ versions.yml ]
}
