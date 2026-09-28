include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record MarkdupBam {
    id:       String
    bam:      Path?
    cram:     Path?
    bai:      Path
    metrics:  Path
    samtools: SamtoolsStatsFiles?
}
