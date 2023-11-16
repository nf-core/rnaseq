//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    ch_bam
        | SAMTOOLS_SORT
        | { out -> SAMTOOLS_INDEX ( out.bam ) }

    SAMTOOLS_SORT.out.bam
        | join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        | join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        | map { meta, bam, bai, csi ->
            if (bai) {
                [ meta, bam, bai ]
            } else {
                [ meta, bam, csi ]
            }
        }
        | { ch_bam_bai -> BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta ) }

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
