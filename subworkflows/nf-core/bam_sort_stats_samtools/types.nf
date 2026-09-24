include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record SortedBam {
    id:       String
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}
