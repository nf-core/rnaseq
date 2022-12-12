//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {

    take:
    bam   // channel: [ val(meta), [ bam ] ]
    fasta // channel: [ fasta ]
    fai   // channel: [ fai ]

    main:

    ch_versions = Channel.empty()

    PICARD_MARKDUPLICATES ( bam, fasta, fai )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    PICARD_MARKDUPLICATES.out.bam
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
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam     // channel: [ val(meta), [ bam ] ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    stats    = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    flagstat = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    idxstats = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
