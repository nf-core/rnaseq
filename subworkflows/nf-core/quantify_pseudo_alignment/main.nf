//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { PseudoQuantSample                   } from './types'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: [ val(meta2), /path/to/index/ ]
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner            //     val: kallisto or salmon
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode
    skip_merge                //    bool: skip cross-sample merging, run tximport per-sample

    main:

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT (
            reads,
            index.combine(gtf).combine(transcript_fasta).first()
        )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results

        // Salmon writes its log inside the quant directory rather than as a
        // discrete file, so log is null. meta_info.json is an optional output.
        ch_sample_fields = SALMON_QUANT.out.results.map { meta, dir -> [meta.id, meta, dir] }
            .join(SALMON_QUANT.out.json_info.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .map { id, meta, dir, json_info -> [id, [meta: meta, quant_dir: dir, json_info: json_info, log: null]] }
    } else {
        KALLISTO_QUANT (
            reads,
            index.combine(gtf.map { g -> [ g, [] ] }).first(),
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log

        ch_sample_fields = KALLISTO_QUANT.out.results.map { meta, dir -> [meta.id, meta, dir] }
            .join(KALLISTO_QUANT.out.json_info.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .join(KALLISTO_QUANT.out.log.map { meta, f -> [meta.id, f] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
            .map { id, meta, dir, json_info, log -> [id, [meta: meta, quant_dir: dir, json_info: json_info, log: log]] }
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT (
        samplesheet,
        ch_pseudo_results,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner,
        skip_merge
    )

    //
    // Merging yields a single 'all_samples' QuantMerged row, which is broadcast
    // onto every sample; under skip_merge there is one row per sample id.
    //
    ch_sample_results = skip_merge
        ? ch_sample_fields.join(QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.results.map { r -> [r.id, r] }, by: [0], failOnMismatch: true, failOnDuplicate: true)
        : ch_sample_fields.combine(QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.results)

    ch_sample_results = ch_sample_results.map { id, fields, quant_merged ->
        record(
            id:           id,
            meta:         fields.meta,
            quant_dir:    fields.quant_dir,
            json_info:    fields.json_info,
            log:          fields.log,
            quant_merged: quant_merged
        )
    }

    emit:
    results                       = ch_pseudo_results                                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                                              // channel: [ val(meta), files_for_multiqc ]
    tx2gene                       = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tx2gene                // channel: [ val(meta), tx2gene.tsv ]
    tx2gene_augmented             = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tx2gene_augmented      // channel: [ val(meta), tx2gene_augmented.tsv ]

    tpm_gene                      = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_gene               //    path: *gene_tpm.tsv
    counts_gene                   = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene            //    path: *gene_counts.tsv
    lengths_gene                  = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_gene           //    path: *gene_lengths.tsv
    counts_gene_length_scaled     = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled            = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_gene_scaled     //    path: *gene_counts_scaled.tsv
    tpm_transcript                = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.tpm_transcript         //    path: *transcript_tpm.tsv
    counts_transcript             = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.counts_transcript      //    path: *transcript_counts.tsv
    lengths_transcript            = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.lengths_transcript     //    path: *transcript_lengths.tsv

    merged_gene_rds_unified       = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_gene_rds       //    path: *.rds
    merged_transcript_rds_unified = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out.merged_transcript_rds //    path: *.rds

    sample_results                = ch_sample_results                                              // channel: PseudoQuantSample
}
