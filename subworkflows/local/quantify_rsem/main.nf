//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION } from '../../../modules/nf-core/rsem/calculateexpression'
include { RSEM_MERGE_COUNTS        } from '../../../modules/local/rsem_merge_counts'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../nf-core/bam_sort_stats_samtools'
include { SENTIEON_RSEMCALCULATEEXPRESSION } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

workflow QUANTIFY_RSEM {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    index        // channel: /path/to/rsem/index/
    fasta        // channel: [ val(meta), [ fasta ] ]
    use_sentieon_star // boolean: determines whether RSEM is run with Sentieon accelerated STAR

    main:

    ch_versions = Channel.empty()

    //
    // Quantify reads with RSEM
    //
    ch_rsem_out = Channel.empty
    if (use_sentieon_star){
        SENTIEON_RSEMCALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = SENTIEON_RSEMCALCULATEEXPRESSION
    } else {
        RSEM_CALCULATEEXPRESSION ( reads, index )
        ch_rsem_out = RSEM_CALCULATEEXPRESSION
    }

    ch_bam_star = ch_rsem_out.out.bam_star
    ch_counts_gene = ch_rsem_out.out.counts_gene
    ch_counts_transcript = ch_rsem_out.out.counts_transcript
    ch_stat = ch_rsem_out.out.stat
    ch_logs = ch_rsem_out.out.logs
    ch_bam_genome = ch_rsem_out.out.bam_genome
    ch_bam_transcript = ch_rsem_out.out.bam_transcript
    ch_versions = ch_versions.mix(ch_rsem_out.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( ch_bam_star, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // Merge counts across samples
    //
    RSEM_MERGE_COUNTS (
        ch_counts_gene.collect{it[1]},       // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        ch_counts_transcript.collect{it[1]}
    )
    ch_versions = ch_versions.mix(RSEM_MERGE_COUNTS.out.versions)

    emit:
    counts_gene              = ch_counts_gene       // channel: [ val(meta), counts ]
    counts_transcript        = ch_counts_transcript // channel: [ val(meta), counts ]
    stat                     = ch_stat              // channel: [ val(meta), stat ]
    logs                     = ch_logs              // channel: [ val(meta), logs ]
    bam_star                 = ch_bam_star          // channel: [ val(meta), bam ]
    bam_genome               = ch_bam_genome        // channel: [ val(meta), bam ]
    bam_transcript           = ch_bam_transcript    // channel: [ val(meta), bam ]

    bam                      = BAM_SORT_STATS_SAMTOOLS.out.bam                // channel: [ val(meta), [ bam ] ]
    bai                      = BAM_SORT_STATS_SAMTOOLS.out.bai                // channel: [ val(meta), [ bai ] ]
    csi                      = BAM_SORT_STATS_SAMTOOLS.out.csi                // channel: [ val(meta), [ csi ] ]
    stats                    = BAM_SORT_STATS_SAMTOOLS.out.stats              // channel: [ val(meta), [ stats ] ]
    flagstat                 = BAM_SORT_STATS_SAMTOOLS.out.flagstat           // channel: [ val(meta), [ flagstat ] ]
    idxstats                 = BAM_SORT_STATS_SAMTOOLS.out.idxstats           // channel: [ val(meta), [ idxstats ] ]

    merged_counts_gene       = RSEM_MERGE_COUNTS.out.counts_gene              //    path: *.gene_counts.tsv
    merged_tpm_gene          = RSEM_MERGE_COUNTS.out.tpm_gene                 //    path: *.gene_tpm.tsv
    merged_counts_transcript = RSEM_MERGE_COUNTS.out.counts_transcript        //    path: *.transcript_counts.tsv
    merged_tpm_transcript    = RSEM_MERGE_COUNTS.out.tpm_transcript           //    path: *.transcript_tpm.tsv

    versions                 = ch_versions                                    // channel: [ versions.yml ]
}
