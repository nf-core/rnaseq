include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/hisat2/index
    splicesites // channel: /path/to/genome.splicesites.txt
    ch_fasta    // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()


    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN.config.ext.args   = '--met-stderr --new-summary --dta'
    HISAT2_ALIGN.config.publishDir = [
        [
            path: "${params.outdir}/${params.aligner}/log",
            pattern: '*.log'
        ],
        [
            path: "${params.outdir}/${params.aligner}",
            pattern: '*.bam',
            enabled: params.save_align_intermeds
        ],
        [
            path: "${params.outdir}/${params.aligner}/unmapped",
            pattern: '*.fastq.gz',
            enabled: params.save_unaligned
        ]
    ]
    HISAT2_ALIGN ( reads, index, splicesites )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    sort_ext_prefix = { "${meta.id}.sorted" }
    sort_publish_dir = [
        path: "${params.outdir}/${params.aligner}",
        pattern: "*.bam",
        enabled: ( !params.with_umi && params.skip_markduplicates ) ||
            params.save_align_intermeds ||
            params.skip_markduplicates
    ]
    index_ext_args = params.bam_csi_index ? '-c' : ''
    index_publish_dir = [
        path: "${params.outdir}/${params.aligner}",
        pattern: "*.{bai,csi}",
        enabled: ( !params.with_umi && params.skip_markduplicates ) ||
            params.save_align_intermeds ||
            params.skip_markduplicates
    ]
    stats_ext_prefix = { "${meta.id}.sorted.bam" }
    stats_publish_dir = [
        path: "${params.outdir}/${params.aligner}/samtools_stats",
        pattern: "*.{stats,flagstat,idxstats}"
    ]
    BAM_SORT_STATS_SAMTOOLS (
        HISAT2_ALIGN.out.bam,
        ch_fasta,
        sort_ext_prefix,
        sort_publish_dir,
        index_ext_args,
        index_publish_dir,
        stats_ext_prefix,
        stats_publish_dir
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)


    emit:
    orig_bam = HISAT2_ALIGN.out.bam                 // channel: [ val(meta), bam   ]
    summary  = HISAT2_ALIGN.out.summary             // channel: [ val(meta), log   ]
    fastq    = HISAT2_ALIGN.out.fastq               // channel: [ val(meta), fastq ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}

