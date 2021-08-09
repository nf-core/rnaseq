//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

params.sort_options  = [:]
params.index_options = [:]
params.stats_options = [:]

include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'  addParams( options: params.sort_options  )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: params.index_options )
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools'                              addParams( options: params.stats_options )

workflow BAM_SORT_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
    SAMTOOLS_SORT      ( ch_bam )
    SAMTOOLS_INDEX     ( SAMTOOLS_SORT.out.bam )

    SAMTOOLS_SORT.out.bam
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

    BAM_STATS_SAMTOOLS ( ch_bam_bai )

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    version  = SAMTOOLS_SORT.out.version       //    path: versions.yml
}
