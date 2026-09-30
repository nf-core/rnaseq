//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SamtoolsStats     } from './types'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta) ]

    main:
    SAMTOOLS_STATS(ch_bam_bai, ch_fasta_fai)

    SAMTOOLS_FLAGSTAT(ch_bam_bai)

    SAMTOOLS_IDXSTATS(ch_bam_bai)

    ch_results = SAMTOOLS_STATS.out
        .join(SAMTOOLS_FLAGSTAT.out, by: 'meta')
        .join(SAMTOOLS_IDXSTATS.out, by: 'meta')
        .map { r ->
            record(id: r.meta.id, stats: r.stats, flagstat: r.flagstat, idxstats: r.idxstats)
        }

    emit:
    stats    = SAMTOOLS_STATS.out.map { r -> [r.meta, r.stats] } // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.map { r -> [r.meta, r.flagstat] } // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.map { r -> [r.meta, r.idxstats] } // channel: [ val(meta), path(idxstats) ]
    results  = ch_results // channel: SamtoolsStats
}
