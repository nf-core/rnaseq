//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'

record RsemMerge {
    counts_gene:       Path
    tpm_gene:          Path
    counts_transcript: Path
    tpm_transcript:    Path
    genes_long:        Path
    isoforms_long:     Path
}

// The QuantMerged record emitted by QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT plus
// rsem_merge. nf-core tools treats every name included from a sibling
// subworkflow as a subworkflow dependency, so the fields are redeclared.
record RsemQuantMerged {
    id:                        String
    meta:                      Map
    tpm_gene:                  Path
    counts_gene:               Path
    lengths_gene:              Path
    counts_gene_length_scaled: Path
    counts_gene_scaled:        Path
    tpm_transcript:            Path
    counts_transcript:         Path
    lengths_transcript:        Path
    tx2gene:                   Path
    tx2gene_augmented:         Path
    merged_gene_rds:           Path?
    merged_transcript_rds:     Path?
    rsem_merge:                RsemMerge?
}

record RsemQuantSample {
    id:                String
    meta:              Map
    counts_gene:       Path
    counts_transcript: Path
    stat:              Path
    log:               Path?
    bam_star:          Path?
    bam_genome:        Path?
    bam_transcript:    Path?
    quant_merged:      RsemQuantMerged
}

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

    // The log and BAMs are optional module outputs, the BAMs depending on ext.args.
    ch_sample_fields = ch_counts_gene.map { meta, f -> [meta.id, meta, f] }
        .join(ch_counts_transcript.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
        .join(ch_stat.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
        .join(ch_logs.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
        .join(ch_rsem_out.out.bam_star.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
        .join(ch_rsem_out.out.bam_genome.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
        .join(ch_rsem_out.out.bam_transcript.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
        .map { id, meta, counts_gene, counts_transcript, stat, log, bam_star, bam_genome, bam_transcript ->
            [id, [
                meta:              meta,
                counts_gene:       counts_gene,
                counts_transcript: counts_transcript,
                stat:              stat,
                log:               log,
                bam_star:          bam_star,
                bam_genome:        bam_genome,
                bam_transcript:    bam_transcript
            ]]
        }

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
        ch_merged_counts_gene       = CUSTOM_RSEMMERGECOUNTS.out.counts_gene
        ch_merged_tpm_gene          = CUSTOM_RSEMMERGECOUNTS.out.tpm_gene
        ch_merged_counts_transcript = CUSTOM_RSEMMERGECOUNTS.out.counts_transcript
        ch_merged_tpm_transcript    = CUSTOM_RSEMMERGECOUNTS.out.tpm_transcript
        ch_merged_genes_long        = CUSTOM_RSEMMERGECOUNTS.out.genes_long
        ch_merged_isoforms_long     = CUSTOM_RSEMMERGECOUNTS.out.isoforms_long

        ch_rsem_merge = ch_merged_counts_gene.map { meta, f -> [meta.id, f] }
            .join(ch_merged_tpm_gene.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .join(ch_merged_counts_transcript.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .join(ch_merged_tpm_transcript.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .join(ch_merged_genes_long.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .join(ch_merged_isoforms_long.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .map { id, counts_gene, tpm_gene, counts_transcript, tpm_transcript, genes_long, isoforms_long ->
                [id, record(
                    counts_gene:       counts_gene,
                    tpm_gene:          tpm_gene,
                    counts_transcript: counts_transcript,
                    tpm_transcript:    tpm_transcript,
                    genes_long:        genes_long,
                    isoforms_long:     isoforms_long
                )]
            }
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

    //
    // Merging yields a single 'all_samples' row from both tximport and
    // CUSTOM_RSEMMERGECOUNTS, joined on that id and broadcast onto every
    // sample. Under skip_merge tximport runs per sample and there is no
    // RSEM merge.
    //
    if (skip_merge) {
        ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.results
            .map { r -> [r.id, r + record(rsem_merge: null)] }
        ch_results = ch_sample_fields
            .join(ch_quant_merged, by: [0], failOnMismatch: true, failOnDuplicate: true)
    } else {
        ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.results
            .map { r -> [r.id, r] }
            .join(ch_rsem_merge, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .map { _id, r, rsem_merge -> r + record(rsem_merge: rsem_merge) }
        ch_results = ch_sample_fields.combine(ch_quant_merged)
    }

    ch_results = ch_results.map { id, fields, quant_merged ->
        record(
            id:                id,
            meta:              fields.meta,
            counts_gene:       fields.counts_gene,
            counts_transcript: fields.counts_transcript,
            stat:              fields.stat,
            log:               fields.log,
            bam_star:          fields.bam_star,
            bam_genome:        fields.bam_genome,
            bam_transcript:    fields.bam_transcript,
            quant_merged:      quant_merged
        )
    }

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
    results                   = ch_results                                                           // channel: RsemQuantSample
}
