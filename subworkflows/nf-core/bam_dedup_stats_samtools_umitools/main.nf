//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP     } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai/csi) ]
    val_get_dedup_stats // boolean: true/false
    dedup_ext_prefix
    index_ext_args
    index_ext_prefix
    index_publish_dir
    stats_ext_prefix

    main:

    ch_versions = Channel.empty()

    //
    // UMI-tools dedup
    //
    UMITOOLS_DEDUP.config.ext.args   = { [
        meta.single_end                 ? '' : '--unpaired-reads=discard --chimeric-pairs=discard',
        params.umitools_grouping_method ? "--method='${params.umitools_grouping_method}'" : '',
        params.umitools_umi_separator   ? "--umi-separator='${params.umitools_umi_separator}'" : ''
    ].join(' ').trim() }
    UMITOOLS_DEDUP.config.ext.prefix = dedup_ext_prefix
    UMITOOLS_DEDUP.config.publishDir = [
        [
            path: "${params.outdir}/${params.aligner}/umitools",
            pattern: '*.tsv'
        ],
        [
            path: "${params.outdir}/${params.aligner}",
            pattern: '*.bam',
            enabled: (
                params.save_align_intermeds ||
                params.save_umi_intermeds
            )
        ]
    ]
    UMITOOLS_DEDUP ( ch_bam_bai, val_get_dedup_stats )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX.config.ext.args   = index_ext_args
    SAMTOOLS_INDEX.config.ext.prefix = index_ext_prefix
    SAMTOOLS_INDEX.config.publishDir = index_publish_dir
    SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_bai_dedup = UMITOOLS_DEDUP.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }

    BAM_STATS_SAMTOOLS.config.ext.prefix = stats_ext_prefix
    BAM_STATS_SAMTOOLS.config.publishDir = [
        path: "${params.outdir}/${params.aligner}/samtools_stats",
        pattern: '*.{stats,flagstat,idxstats}'
    ]
    BAM_STATS_SAMTOOLS ( ch_bam_bai_dedup, [ [:], [] ] )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = UMITOOLS_DEDUP.out.bam          // channel: [ val(meta), path(bam) ]

    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), path(csi) ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                     // channel: [ path(versions.yml) ]
}
