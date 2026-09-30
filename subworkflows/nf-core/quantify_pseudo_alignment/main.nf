nextflow.enable.types = true

//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { Reads                               } from './types'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<Reads>
    index: Value<Tuple<Map, Path>>
    transcript_fasta: Value<Path>
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner: String // kallisto or salmon
    kallisto_quant_fraglen: Integer? // Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd: Integer? // Estimated standard error for fragment length required by Kallisto in single-end mode
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample

    main:

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        ch_salmon = SALMON_QUANT(ch_samples, index.combine(gtf).combine(transcript_fasta))

        // Salmon writes its log inside the quant directory rather than as a
        // discrete file, so log is null. meta_info.json is an optional output.
        ch_sample_results = ch_salmon.map { r ->
            record(id: r.id, meta: r.meta, quant_dir: r.quant_dir, json_info: r.json_info, log: null, multiqc: r.quant_dir)
        }
    } else {
        ch_kallisto = KALLISTO_QUANT(
            ch_samples,
            index.combine(gtf).map { meta, idx, gtf_file -> tuple(meta, idx, gtf_file, null) },
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_sample_results = ch_kallisto.map { r ->
            record(id: r.id, meta: r.meta, quant_dir: r.quant_dir, json_info: r.json_info, log: r.log, multiqc: r.log)
        }
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_sample_results.map { r -> record(id: r.id, meta: r.meta, quants: [ r.quant_dir ]) },
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner,
        skip_merge
    )

    emit:
    samples = ch_sample_results // per sample
    merged  = ch_quant_merged   // one row per sample under skip_merge, a single 'all_samples' row otherwise
}
