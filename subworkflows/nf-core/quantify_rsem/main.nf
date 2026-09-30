nextflow.enable.types = true

//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { Reads; RsemQuantified               } from './types'

workflow QUANTIFY_RSEM {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<Reads> // FASTQ or BAM files
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
    if (skip_merge) {
        ch_rsem_merge = ch_rsem.map { r -> record(id: r.id, rsem_merge: null) }
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
    ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_rsem.map { r -> record(id: r.id, meta: r.meta, quants: [ r.counts_transcript ]) },
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        'rsem',
        skip_merge
    )

    //
    // Per-sample rows carry the RSEM outputs; merged rows carry the tximport (and, when
    // merging, RSEM merge) outputs. A merged row has the id of its tximport run, which is a
    // sample id under skip_merge and 'all_samples' otherwise.
    //
    ch_sample_rows = ch_rsem.map { r -> record(id: r.id, meta: r.meta, sample: r, merged: null, rsem_merge: null) }
    ch_merged_rows = ch_quant_merged
        .map { m -> record(id: m.id, meta: m.meta, merged: m) }
        .join(ch_rsem_merge, by: 'id')
        .map { r -> r + record(sample: null) }

    ch_results = ch_sample_rows.mix(ch_merged_rows)

    emit:
    ch_results
}
