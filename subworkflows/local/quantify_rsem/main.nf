//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION } from '../../../modules/nf-core/rsem/calculateexpression'
include { RSEM_MERGE_COUNTS        } from '../../../modules/local/rsem_merge_counts'
include { SENTIEON_RSEMCALCULATEEXPRESSION } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

workflow QUANTIFY_RSEM {
    take:
    reads             // channel: [ val(meta), [ reads ] ] - FASTQ or BAM files
    index             // channel: /path/to/rsem/index/
    use_sentieon_star // boolean: determines whether RSEM is run with Sentieon accelerated STAR

    main:

    ch_versions = Channel.empty()

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

    ch_counts_gene = ch_rsem_out.out.counts_gene
    ch_counts_transcript = ch_rsem_out.out.counts_transcript
    ch_stat = ch_rsem_out.out.stat
    ch_logs = ch_rsem_out.out.logs
    ch_versions = ch_versions.mix(ch_rsem_out.out.versions.first())

    //
    // Merge counts across samples
    //
    RSEM_MERGE_COUNTS (
        ch_counts_gene.collect{it[1]},       // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        ch_counts_transcript.collect{it[1]}
    )
    ch_versions = ch_versions.mix(RSEM_MERGE_COUNTS.out.versions)

    emit:
    counts_gene              = ch_counts_gene                                 // channel: [ val(meta), counts ]
    counts_transcript        = ch_counts_transcript                           // channel: [ val(meta), counts ]
    stat                     = ch_stat                                        // channel: [ val(meta), stat ]
    logs                     = ch_logs                                        // channel: [ val(meta), logs ]

    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.tsv
    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.tsv
    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.tsv
    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.tsv

    versions                 = ch_versions                                    // channel: [ versions.yml ]
}
