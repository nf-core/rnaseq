//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {

    take:
    ch_bam   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ path(fasta) ]
    ch_fai   // channel: [ path(fai) ]

    main:

    ch_versions = Channel.empty()

    PICARD_MARKDUPLICATES ( ch_bam, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_bai = PICARD_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map{meta, bam, bai, csi ->
            if (bai) [ meta, bam, bai ]
            else [ meta, bam, csi ]
        }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), path(bam) ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), path(csi) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}

workflow {
    main:
    def ch_bam = Channel.empty()
    def ch_fasta = Channel.empty()
    def ch_fai = Channel.empty()

    BAM_MARKDUPLICATES_PICARD( ch_bam, ch_fasta, ch_fai )

    publish:
    BAM_MARKDUPLICATES_PICARD.out.bam       >> 'picard'
    BAM_MARKDUPLICATES_PICARD.out.metrics   >> 'picard/metrics'
    BAM_MARKDUPLICATES_PICARD.out.bai       >> 'picard'
    BAM_MARKDUPLICATES_PICARD.out.csi       >> 'picard'
    BAM_MARKDUPLICATES_PICARD.out.stats     >> 'picard/samtools_stats'
    BAM_MARKDUPLICATES_PICARD.out.flagstat  >> 'picard/samtools_stats'
    BAM_MARKDUPLICATES_PICARD.out.idxstats  >> 'picard/samtools_stats'
}
