nextflow.enable.types = true

//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { ReadsInput; SalmonQuantSample; KallistoQuantSample } from '../../../modules/nf-core/types'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<ReadsInput>
    index: Value<Path?>
    transcript_fasta: Value<Path?>
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
    def ch_salmon: Channel<SalmonQuantSample>     = channel.empty()
    def ch_kallisto: Channel<KallistoQuantSample> = channel.empty()
    if (pseudo_aligner == 'salmon') {
        ch_salmon = SALMON_QUANT(ch_samples, index, gtf, transcript_fasta)
        ch_quant_dirs = ch_salmon.map { r -> record(id: r.id, meta: r.meta, quants: [ r.quant_dir ]) }
    } else {
        ch_kallisto = KALLISTO_QUANT(
            ch_samples,
            index,
            gtf,
            null,
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_quant_dirs = ch_kallisto.map { r -> record(id: r.id, meta: r.meta, quants: [ r.quant_dir ]) }
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_quant_dirs,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner,
        skip_merge
    )

    emit:
    salmon   = ch_salmon        // per sample, when pseudo_aligner is salmon
    kallisto = ch_kallisto      // per sample, when pseudo_aligner is kallisto
    merged   = ch_quant_merged  // one row per sample under skip_merge, a single 'all_samples' row otherwise
}
