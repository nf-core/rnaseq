// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record Hisat2Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

record Hisat2Logs {
    summary: Path
}

record Hisat2Aligned {
    id:       String
    meta:     Map
    aligner:  String
    raw_bams: List<Path>
    unmapped: List<Path>
    hisat2:   Hisat2Logs
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}
