nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { Bam                } from '../../local/types'
include { SamtoolsStats      } from './types'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<Bam>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    ch_stats = SAMTOOLS_STATS(ch_bam_bai, ch_fasta, ch_fai)
    ch_flagstat = SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_idxstats = SAMTOOLS_IDXSTATS(ch_bam_bai)

    ch_stats_flagstat = ch_stats.join(ch_flagstat, by: 'id', remainder: true)
    ch_stats_flagstat.subscribe { r ->
        if( r.stats == null || r.flagstat == null ) {
            error "Sample '${r.id}' is missing its samtools flagstat result"
        }
    }
    ch_stats_idxstats = ch_stats_flagstat
        .filter { r -> r.stats != null && r.flagstat != null }
        .join(ch_idxstats, by: 'id', remainder: true)
    ch_stats_idxstats.subscribe { r ->
        if( r.flagstat == null || r.idxstats == null ) {
            error "Sample '${r.id}' is missing its samtools idxstats result"
        }
    }
    ch_results = ch_stats_idxstats
        .filter { r -> r.flagstat != null && r.idxstats != null }
        .map { r -> record(id: r.id, meta: r.meta, samtools: record(stats: r.stats, flagstat: r.flagstat, idxstats: r.idxstats)) }

    emit:
    ch_results
}
