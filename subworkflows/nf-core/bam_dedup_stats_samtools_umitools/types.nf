include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record UmitoolsDedupStats {
    edit_distance:    Path
    per_umi:          Path
    umi_per_position: Path
}

record UmitoolsDedupBam {
    id:        String
    bam:       Path
    bai:       Path
    dedup_log: Path
    samtools:  SamtoolsStatsFiles
    tsv:       UmitoolsDedupStats?
}
