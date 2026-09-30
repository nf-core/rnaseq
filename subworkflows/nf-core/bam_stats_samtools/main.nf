nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { Bam; SamtoolsStats; SamtoolsStatsResult; SamtoolsFlagstatResult; SamtoolsIdxstatsResult } from '../../../modules/nf-core/types'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<Bam>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    def ch_stats: Channel<SamtoolsStatsResult> = SAMTOOLS_STATS(ch_bam_bai, ch_fasta, ch_fai)
    def ch_flagstat: Channel<SamtoolsFlagstatResult> = SAMTOOLS_FLAGSTAT(ch_bam_bai)
    def ch_idxstats: Channel<SamtoolsIdxstatsResult> = SAMTOOLS_IDXSTATS(ch_bam_bai)

    def ch_results: Channel<SamtoolsStats> = ch_stats
        .join(ch_flagstat, by: 'id')
        .join(ch_idxstats, by: 'id')
        .map { r -> record(id: r.id, meta: r.meta, samtools: record(stats: r.stats, flagstat: r.flagstat, idxstats: r.idxstats)) }

    emit:
    ch_results
}
