nextflow.enable.types = true

//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { ReadsInput; QuantMerged } from '../../../modules/nf-core/types'
include { RsemMergeSample } from '../../../modules/nf-core/custom/rsemmergecounts/main'
include { RsemQuantSample } from '../../../modules/nf-core/rsem/calculateexpression/main'

workflow QUANTIFY_RSEM {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<ReadsInput> // FASTQ or BAM files
    index: Value<Path> // RSEM index
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    use_sentieon_star: Boolean // use Sentieon-accelerated STAR (FASTQ mode only)
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample

    main:

    //
    // Quantify reads with RSEM
    //
    def ch_rsem: Channel<RsemQuantSample>
    if (use_sentieon_star) {
        ch_rsem = SENTIEON_RSEMCALCULATEEXPRESSION(ch_samples, index)
    } else {
        ch_rsem = RSEM_CALCULATEEXPRESSION(ch_samples, index)
    }

    //
    // Merge counts across samples (only when not skipping merge)
    //
    // Sorted by name for a stable cache key; the script globs the
    // staged directory, so order doesn't affect output. The emptiness
    // check skips CUSTOM_RSEMMERGECOUNTS when there are no samples, since
    // collect emits [] rather than nothing on an empty channel.
    //
    def ch_rsem_merge: Channel<RsemMergeSample>
    if (skip_merge) {
        ch_rsem_merge = channel.empty()
    } else {
        ch_rsem_merge = CUSTOM_RSEMMERGECOUNTS(
            ch_rsem
                .collect()
                .flatMap { rs ->
                    rs.isEmpty()
                        ? []
                        : [ record(
                            id: 'all_samples',
                            meta: [id: 'all_samples'],
                            genes: rs.collect { r -> r.counts_gene }.toSorted { f -> f.name },
                            isoforms: rs.collect { r -> r.counts_transcript }.toSorted { f -> f.name }
                        ) ]
                }
        ).map { r -> record(id: r.id, rsem_merge: r.rsem_merge) }
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    def ch_quant_merged: Channel<QuantMerged> = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_rsem.map { r -> record(id: r.id, meta: r.meta, quants: [ r.counts_transcript ]) },
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        'rsem',
        skip_merge,
        null
    )

    emit:
    samples:    Channel<RsemQuantSample> = ch_rsem        // per sample
    merged:     Channel<QuantMerged>     = ch_quant_merged // one row per sample under skip_merge, a single 'all_samples' row otherwise
    rsem_merge: Channel<RsemMergeSample> = ch_rsem_merge  // a single 'all_samples' row; empty under skip_merge
}
