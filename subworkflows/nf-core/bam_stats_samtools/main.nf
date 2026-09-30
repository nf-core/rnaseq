nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { BamBaiInput; SamtoolsStats } from '../../../modules/nf-core/types'
include { SamtoolsFlagstatResult } from '../../../modules/nf-core/samtools/flagstat/main'
include { SamtoolsIdxstatsResult } from '../../../modules/nf-core/samtools/idxstats/main'
include { SamtoolsStatsResult } from '../../../modules/nf-core/samtools/stats/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<BamBaiInput>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    def ch_stats: Channel<SamtoolsStatsResult> = SAMTOOLS_STATS(ch_bam_bai, ch_fasta, ch_fai)
    def ch_flagstat: Channel<SamtoolsFlagstatResult> = SAMTOOLS_FLAGSTAT(ch_bam_bai)
    def ch_idxstats: Channel<SamtoolsIdxstatsResult> = SAMTOOLS_IDXSTATS(ch_bam_bai)

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
    def ch_results: Channel<SamtoolsStats> = ch_stats_idxstats
        .filter { r -> r.flagstat != null && r.idxstats != null }
        .map { r -> record(id: r.id, meta: r.meta, samtools: record(stats: r.stats, flagstat: r.flagstat, idxstats: r.idxstats)) }

    emit:
    ch_results
}
