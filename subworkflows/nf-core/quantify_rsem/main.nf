//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { RsemQuantSample                     } from './types'

workflow QUANTIFY_RSEM {
    take:
    samplesheet           // channel: [ val(meta), /path/to/samplesheet ]
    reads                 // channel: [ val(meta), [ reads ] ] - FASTQ or BAM files
    index                 // channel: /path/to/rsem/index/
    gtf                   // channel: /path/to/genome.gtf
    gtf_id_attribute      //     val: GTF gene ID attribute
    gtf_extra_attribute   //     val: GTF alternative gene attribute (e.g. gene_name)
    use_sentieon_star     // boolean: use Sentieon-accelerated STAR (FASTQ mode only)
    skip_merge            //    bool: skip cross-sample merging, run tximport per-sample

    main:

    //
    // Quantify reads with RSEM
    //
    ch_rsem = null
    if (use_sentieon_star) {
        SENTIEON_RSEMCALCULATEEXPRESSION ( reads, index )
        ch_rsem = SENTIEON_RSEMCALCULATEEXPRESSION.out
    } else {
        RSEM_CALCULATEEXPRESSION ( reads, index )
        ch_rsem = RSEM_CALCULATEEXPRESSION.out
    }

    // The log and BAMs are optional module outputs, the BAMs depending on ext.args.
    ch_counts_gene       = ch_rsem.map { r -> [ r.meta, r.counts_gene ] }
    ch_counts_transcript = ch_rsem.map { r -> [ r.meta, r.counts_transcript ] }
    ch_stat              = ch_rsem.map { r -> [ r.meta, r.stat ] }
    ch_logs              = ch_rsem.filter { r -> r.log }.map { r -> [ r.meta, r.log ] }

    //
    // Merge counts across samples (only when not skipping merge)
    //
    ch_merged_counts_gene       = channel.empty()
    ch_merged_tpm_gene          = channel.empty()
    ch_merged_counts_transcript = channel.empty()
    ch_merged_tpm_transcript    = channel.empty()
    ch_merged_genes_long        = channel.empty()
    ch_merged_isoforms_long     = channel.empty()
    ch_rsem_merge               = channel.empty()

    if (!skip_merge) {
        //
        // Sorted by name for a stable cache key; the script globs the
        // staged directory, so order doesn't affect output. Filtered to
        // skip CUSTOM_RSEMMERGECOUNTS when there are no samples, since
        // toSortedList() emits [] rather than nothing on an empty channel.
        //
        CUSTOM_RSEMMERGECOUNTS (
            ch_counts_gene
                .toSortedList { a, b -> a[1].name <=> b[1].name }
                .filter { sorted -> sorted.size() > 0 }
                .map { sorted -> [ ['id': 'all_samples'], sorted.collect { it[1] } ] },
            ch_counts_transcript
                .toSortedList { a, b -> a[1].name <=> b[1].name }
                .filter { sorted -> sorted.size() > 0 }
                .map { sorted -> sorted.collect { it[1] } }
        )
        ch_merged_counts_gene       = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.counts_gene ] }
        ch_merged_tpm_gene          = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.tpm_gene ] }
        ch_merged_counts_transcript = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.counts_transcript ] }
        ch_merged_tpm_transcript    = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.tpm_transcript ] }
        ch_merged_genes_long        = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.genes_long ] }
        ch_merged_isoforms_long     = CUSTOM_RSEMMERGECOUNTS.out.map { r -> [ [id: r.id], r.rsem_merge.isoforms_long ] }

        ch_rsem_merge = CUSTOM_RSEMMERGECOUNTS.out
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT (
        samplesheet,
        ch_counts_transcript,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        'rsem',
        skip_merge
    )

    emit:
    // Per-sample outputs
    raw_counts_gene          = ch_counts_gene                                                        // channel: [ val(meta), counts ]
    raw_counts_transcript    = ch_counts_transcript                                                  // channel: [ val(meta), counts ]
    stat                     = ch_stat                                                               // channel: [ val(meta), stat ]
    logs                     = ch_logs                                                               // channel: [ val(meta), logs ]

    // RSEM merge outputs
    merged_counts_gene       = ch_merged_counts_gene                                                 // channel: [ val(meta), counts ]
    merged_tpm_gene          = ch_merged_tpm_gene                                                    // channel: [ val(meta), tpm ]
    merged_counts_transcript = ch_merged_counts_transcript                                           // channel: [ val(meta), counts ]
    merged_tpm_transcript    = ch_merged_tpm_transcript                                              // channel: [ val(meta), tpm ]
    merged_genes_long        = ch_merged_genes_long                                                  // channel: [ val(meta), genes_long ]
    merged_isoforms_long     = ch_merged_isoforms_long                                               // channel: [ val(meta), isoforms_long ]

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

    // Per-sample record
    results                   = ch_rsem                                                              // channel: RsemQuantSample

    // Merged quantification: a single row when merging, one per sample under skip_merge
    quant_merged              = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.results                      // channel: QuantMerged

    // RSEM merge counts, keyed by id; empty channel under skip_merge
    rsem_merge                = ch_rsem_merge                                                        // channel: record(id, rsem_merge: RsemMerge)
}
