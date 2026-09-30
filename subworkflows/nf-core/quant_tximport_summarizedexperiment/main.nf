nextflow.enable.types = true

//
// Quantification post-processing with tximport and SummarizedExperiment
//

include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_UNIFIED       } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT_UNIFIED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { QuantsInput; QuantMerged } from '../../../modules/nf-core/types'
include { CustomTx2geneResult } from '../../../modules/nf-core/custom/tx2gene/main'
include { TximetaTximportResult } from '../../../modules/nf-core/tximeta/tximport/main'

workflow QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT {
    take:
    samplesheet: Value<Path>
    ch_quants: Channel<QuantsInput> // per-sample quantification results
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    quant_type: String // 'salmon', 'kallisto', or 'rsem'
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample

    main:

    //
    // Create tx2gene mapping from GTF + a single sample's quantification files.
    // tx2gene only uses quant files to discover which GTF attribute corresponds
    // to transcript IDs, so a single sample suffices. This assumes all samples
    // were quantified against the same transcriptome, avoiding the need to run
    // tx2gene independently per sample. If a use case arises requiring mixed
    // transcriptomes, this will need to be revisited.
    //
    // The sample is selected by sorting on the staged results name so the same
    // one is chosen on every run, keeping the CUSTOM_TX2GENE cache key stable
    // across -resume. Picking by arrival order would vary between runs and
    // invalidate the cache for tx2gene and everything downstream of it. The
    // empty-list guard keeps tx2gene from running when no quant results are
    // supplied, since collect still emits an empty list in that case.
    //
    ch_tx2gene_quants = ch_quants
        .collect()
        .flatMap { rs -> rs.isEmpty() ? [] : [ record(id: 'tx2gene', meta: [:], quants: [ rs.collect { r -> r.quants[0] }.toSorted { f -> f.name }.first() ]) ] }

    def ch_tx2gene: Channel<CustomTx2geneResult> = CUSTOM_TX2GENE(
        ch_tx2gene_quants,
        gtf,
        quant_type,
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_tx2gene_file = ch_tx2gene.collect().map { rs -> rs.isEmpty() ? null : rs.toSorted { r -> r.id }.first().tx2gene }

    //
    // Import and summarize quantifications with tximport
    // In per-sample mode, run once per sample instead of collecting all
    //
    // Sorted by name for a stable cache key; the R script derives sample
    // identity from staged file names, not list position. Filtered to
    // skip TXIMETA_TXIMPORT when there are no samples, since
    // collect emits [] rather than nothing on an empty channel.
    //
    if (skip_merge) {
        ch_tximport_input = ch_quants
    } else {
        ch_tximport_input = ch_quants
            .collect()
            .flatMap { rs -> rs.isEmpty() ? [] : [ record(id: 'all_samples', meta: [id: 'all_samples'], quants: rs.collect { r -> r.quants[0] }.toSorted { f -> f.name }) ] }
    }

    def ch_tximport: Channel<TximetaTximportResult> = TXIMETA_TXIMPORT(ch_tximport_input, ch_tx2gene_file, quant_type)

    //
    // Build SummarizedExperiment objects (only when merging)
    //
    if (skip_merge) {
        ch_se_gene       = ch_tximport.map { r -> record(id: r.id, merged_gene_rds: null) }
        ch_se_transcript = ch_tximport.map { r -> record(id: r.id, merged_transcript_rds: null) }
    } else {
        //
        // Build gene-level SummarizedExperiment
        //
        ch_se_gene = SE_GENE_UNIFIED(
            ch_tximport.map { r ->
                record(id: r.id, meta: r.meta, matrix_files: [ r.counts_gene, r.counts_gene_length_scaled, r.counts_gene_scaled, r.lengths_gene, r.tpm_gene ])
            },
            ch_tx2gene_file,
            samplesheet
        ).map { r -> record(id: r.id, merged_gene_rds: r.rds) }

        //
        // Build transcript-level SummarizedExperiment
        //
        ch_se_transcript = SE_TRANSCRIPT_UNIFIED(
            ch_tximport.map { r ->
                record(id: r.id, meta: r.meta, matrix_files: [ r.counts_transcript, r.lengths_transcript, r.tpm_transcript ])
            },
            ch_tximport.collect().map { rs -> rs.isEmpty() ? null : rs.toSorted { r -> r.id }.first().tx2gene_augmented },
            samplesheet
        ).map { r -> record(id: r.id, merged_transcript_rds: r.rds) }
    }

    //
    // One record per TXIMETA_TXIMPORT row: a single 'all_samples' row when
    // merging, or one per sample under skip_merge. The SE outputs only exist
    // when merging and are null otherwise.
    //
    ch_results = ch_tximport
        .join(ch_se_gene, by: 'id')
        .join(ch_se_transcript, by: 'id')
        .combine(tx2gene: ch_tx2gene_file)

    emit:
    ch_results
}
