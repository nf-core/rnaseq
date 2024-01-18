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

    main:

    ( ch_bam_bai & val_get_dedup_stats )
        | UMITOOLS_DEDUP 
        | { out -> SAMTOOLS_INDEX ( out.bam ) }

    UMITOOLS_DEDUP.out.bam
        | join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        | join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        | map { meta, bam, bai, csi ->
            if (bai) {
                [ meta, bam, bai ]
            } else {
                [ meta, bam, csi ]
            }
        }
        | { ch_bam_bai_dedup -> BAM_STATS_SAMTOOLS ( ch_bam_bai_dedup, [ [:], [] ] ) }

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
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
