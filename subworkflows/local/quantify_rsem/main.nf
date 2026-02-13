nextflow.preview.types = true

//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION } from '../../../modules/nf-core/rsem/calculateexpression'
include { RSEM_MERGE_COUNTS        } from '../../../modules/local/rsem_merge_counts'
include { SENTIEON_RSEMCALCULATEEXPRESSION } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

record RsemResult {
    meta:                      Map
    stat:                      Path?
    logs:                      Path?
    counts_gene:               Path?
    counts_transcript:         Path?
    merged_counts_gene:        Path?
    merged_counts_transcript:  Path?
    merged_tpm_gene:           Path?
    merged_tpm_transcript:     Path?
    merged_genes_long:         Path?
    merged_isoforms_long:      Path?
}

workflow QUANTIFY_RSEM {
    take:
    reads             // channel: [ val(meta), [ reads ] ] - FASTQ or BAM files
    index             // channel: /path/to/rsem/index/
    use_sentieon_star // boolean: determines whether RSEM is run with Sentieon accelerated STAR

    main:

    ch_versions = channel.empty()

    //
    // Quantify reads with RSEM
    //
    ch_rsem_out = null
    if (use_sentieon_star){
        SENTIEON_RSEMCALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = SENTIEON_RSEMCALCULATEEXPRESSION
    } else {
        RSEM_CALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = RSEM_CALCULATEEXPRESSION
    }

    // Extract individual fields from the process record for downstream use
    ch_rsem_result = ch_rsem_out.out
    ch_counts_gene       = ch_rsem_result.map { r -> [r.meta, r.counts_gene] }
    ch_counts_transcript = ch_rsem_result.map { r -> [r.meta, r.counts_transcript] }
    ch_stat = ch_rsem_result.map { r -> [r.meta, r.stat] }
    ch_logs = ch_rsem_result.map { r -> [r.meta, r.logs] }

    //
    // Merge counts across samples
    //
    RSEM_MERGE_COUNTS (
        ch_counts_gene.collect{ tuple -> tuple[1] },
        ch_counts_transcript.collect{ tuple -> tuple[1] }
    )

    // Combine per-sample process records with pipeline-wide aggregate outputs.
    // The RsemMergedResult record lets us combine() once instead of 6 separate times.
    ch_merged = RSEM_MERGE_COUNTS.out

    emit:
    result = ch_rsem_result
        .combine(ch_merged)
        .map { calc, merged ->
            record(
                meta: calc.meta,
                stat: calc.stat, logs: calc.logs,
                counts_gene: calc.counts_gene, counts_transcript: calc.counts_transcript,
                merged_counts_gene: merged.counts_gene, merged_counts_transcript: merged.counts_transcript,
                merged_tpm_gene: merged.tpm_gene, merged_tpm_transcript: merged.tpm_transcript,
                merged_genes_long: merged.genes_long, merged_isoforms_long: merged.isoforms_long
            )
        }
    stat               = ch_stat
    merged_counts_gene = ch_merged.map { r -> r.counts_gene }
    versions           = ch_versions
}
