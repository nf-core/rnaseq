nextflow.preview.types = true

//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'
include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_UNIFIED       } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT_UNIFIED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

record QuantResult {
    meta:                      Map
    results:                   Path?
    tx2gene:                   Path?
    counts_gene:               Path?
    counts_gene_length_scaled: Path?
    counts_gene_scaled:        Path?
    counts_transcript:         Path?
    lengths_gene:              Path?
    lengths_transcript:        Path?
    tpm_gene:                  Path?
    tpm_transcript:            Path?
    merged_gene_rds:           Path?
    merged_transcript_rds:     Path?
}

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: /path/to//index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner            //     val: kallisto or salmon
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode

    main:
    ch_versions = channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT (
            reads,
            index,
            gtf,
            transcript_fasta,
            alignment_mode,
            lib_type
        )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
    } else {
        KALLISTO_QUANT (
            reads,
            index,
            gtf,
            [],
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_pseudo_results = KALLISTO_QUANT.out.map { r -> [r.meta, r.results] }
        ch_pseudo_multiqc = KALLISTO_QUANT.out.map { r -> [r.meta, r.log] }
    }

    CUSTOM_TX2GENE (
        gtf.map { gtf_file -> [ [:], gtf_file ] },
        ch_pseudo_results.collect{ meta_results -> meta_results[1] }.map { results -> [ [:], results ] },
        pseudo_aligner,
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT (
        ch_pseudo_results.collect{ meta_results -> meta_results[1] }.map { results -> [ ['id': 'all_samples'], results ] },
        CUSTOM_TX2GENE.out.tx2gene,
        pseudo_aligner
    )

    // Extract fields directly from TximportResult record (replaces 5 join() calls)
    ch_gene_unified = TXIMETA_TXIMPORT.out
        .map { r -> tuple(r.meta, [r.counts_gene, r.counts_gene_length_scaled, r.counts_gene_scaled, r.lengths_gene, r.tpm_gene]) }

    SE_GENE_UNIFIED (
        ch_gene_unified,
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )
    ch_versions = ch_versions.mix(SE_GENE_UNIFIED.out.versions)

    ch_transcript_unified = TXIMETA_TXIMPORT.out
        .map { r -> tuple(r.meta, [r.counts_transcript, r.lengths_transcript, r.tpm_transcript]) }

    SE_TRANSCRIPT_UNIFIED (
        ch_transcript_unified,
        CUSTOM_TX2GENE.out.tx2gene,
        samplesheet
    )
    ch_versions = ch_versions.mix(SE_TRANSCRIPT_UNIFIED.out.versions)

    // Combine per-sample results with pipeline-wide aggregate outputs.
    // The TximportResult record lets us combine() once instead of 8 separate times.
    ch_tximport = TXIMETA_TXIMPORT.out

    emit:
    result = ch_pseudo_results
        .combine(CUSTOM_TX2GENE.out.tx2gene.map { _meta, path -> path })
        .combine(ch_tximport)
        .combine(SE_GENE_UNIFIED.out.rds.map { _meta, path -> path })
        .combine(SE_TRANSCRIPT_UNIFIED.out.rds.map { _meta, path -> path })
        .map { meta, results, tx2gene, txi, rds_g, rds_t ->
            record(
                meta: meta,
                results: results, tx2gene: tx2gene,
                counts_gene: txi.counts_gene, counts_gene_length_scaled: txi.counts_gene_length_scaled,
                counts_gene_scaled: txi.counts_gene_scaled, counts_transcript: txi.counts_transcript,
                lengths_gene: txi.lengths_gene, lengths_transcript: txi.lengths_transcript,
                tpm_gene: txi.tpm_gene, tpm_transcript: txi.tpm_transcript,
                merged_gene_rds: rds_g, merged_transcript_rds: rds_t
            )
        }
    multiqc                   = ch_pseudo_multiqc
    counts_gene_length_scaled = ch_tximport.map { r -> [r.meta, r.counts_gene_length_scaled] }
    versions                  = ch_versions
}
