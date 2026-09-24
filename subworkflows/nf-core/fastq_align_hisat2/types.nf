include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record Hisat2Logs {
    summary: Path
}

record Hisat2Aligned {
    id:       String
    meta:     Map
    aligner:  String
    orig_bam: Path
    unmapped: List<Path>?
    hisat2:   Hisat2Logs
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}
