include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record UmicollapseDedupBam {
    id:          String
    bam:         Path
    bai:         Path
    dedup_stats: Path
    samtools:    SamtoolsStatsFiles
}
