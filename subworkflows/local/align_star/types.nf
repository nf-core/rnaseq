include { SamtoolsStatsFiles } from '../../nf-core/bam_stats_samtools/types'

record StarLogs {
    log_final:    Path
    log_out:      Path
    log_progress: Path
    tab:          List<Path>?
}

record StarAligned {
    id:                String
    meta:              Map
    aligner:           String
    orig_bam:          List<Path>
    transcriptome_bam: Path?
    unmapped:          List<Path>?
    percent_mapped:    Float
    star:              StarLogs
    bam:               Path
    bai:               Path
    samtools:          SamtoolsStatsFiles
}
