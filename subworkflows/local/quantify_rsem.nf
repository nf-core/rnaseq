//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION } from '../../modules/nf-core/rsem/calculateexpression/main'
include { RSEM_MERGE_COUNTS        } from '../../modules/local/rsem_merge_counts'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow QUANTIFY_RSEM {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/rsem/index/
    fasta // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Quantify reads with RSEM
    //
    RSEM_CALCULATEEXPRESSION.config.ext.args   = [
        '--star',
        '--star-output-genome-bam',
        '--star-gzipped-read-file',
        '--estimate-rspd',
        '--seed 1'
    ].join(' ').trim()
    RSEM_CALCULATEEXPRESSION.config.publishDir = [
        [
            path: "${params.outdir}/${params.aligner}",
            pattern: "*.{stat,results}"
        ],
        [
            path: "${params.outdir}/${params.aligner}",
            pattern: "*.bam",
            enabled: params.save_align_intermeds
        ],
        [
            path: "${params.outdir}/${params.aligner}/log",
            pattern: "*.log"
        ]
    ]
    RSEM_CALCULATEEXPRESSION ( reads, index )
    ch_versions = ch_versions.mix(RSEM_CALCULATEEXPRESSION.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    sort_ext_prefix = { "${meta.id}.sorted" }
    sort_publish_dir = [
        path: "${params.outdir}/${params.aligner}",
        pattern: "*.bam",
        enabled: params.save_align_intermeds ||
            params.skip_markduplicates
    ]
    index_ext_args = params.bam_csi_index ? '-c' : ''
    index_publish_dir = [
        path: "${params.outdir}/${params.aligner}",
        pattern: "*.{bai,csi}",
        enabled: params.save_align_intermeds ||
            params.skip_markduplicates
    ]
    stats_ext_prefix = { "${meta.id}.sorted.bam" }
    stats_publish_dir = [
        path: "${params.outdir}/${params.aligner}/samtools_stats",
        pattern: "*.{stats,flagstat,idxstats}"
    ]
    BAM_SORT_STATS_SAMTOOLS (
        RSEM_CALCULATEEXPRESSION.out.bam_star,
        fasta,
        sort_ext_prefix,
        sort_publish_dir,
        index_ext_args,
        index_publish_dir,
        stats_ext_prefix,
        stats_publish_dir
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // Merge counts across samples
    //
    RSEM_MERGE_COUNTS.config.publishDir = [
        path: "${params.outdir}/${params.aligner}",
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
    RSEM_MERGE_COUNTS (
        RSEM_CALCULATEEXPRESSION.out.counts_gene.collect{it[1]},       // [meta, counts]: Collect the second element (counts files) in the channel across all samples
        RSEM_CALCULATEEXPRESSION.out.counts_transcript.collect{it[1]}
    )
    ch_versions = ch_versions.mix(RSEM_MERGE_COUNTS.out.versions)

    emit:
    counts_gene              = RSEM_CALCULATEEXPRESSION.out.counts_gene       // channel: [ val(meta), counts ]
    counts_transcript        = RSEM_CALCULATEEXPRESSION.out.counts_transcript // channel: [ val(meta), counts ]
    stat                     = RSEM_CALCULATEEXPRESSION.out.stat              // channel: [ val(meta), stat ]
    logs                     = RSEM_CALCULATEEXPRESSION.out.logs              // channel: [ val(meta), logs ]
    bam_star                 = RSEM_CALCULATEEXPRESSION.out.bam_star          // channel: [ val(meta), bam ]
    bam_genome               = RSEM_CALCULATEEXPRESSION.out.bam_genome        // channel: [ val(meta), bam ]
    bam_transcript           = RSEM_CALCULATEEXPRESSION.out.bam_transcript    // channel: [ val(meta), bam ]

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
