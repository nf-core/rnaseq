include { SamtoolsStatsFiles } from '../../nf-core/bam_sort_stats_samtools'

record Bowtie2Logs {
    log: Path
}

record Bowtie2Aligned {
    id:             String
    meta:           Map
    aligner:        String
    orig_bam:       Path
    unmapped:       List<Path>?
    percent_mapped: Float
    bowtie2:        Bowtie2Logs
    bam:            Path
    bai:            Path
    samtools:       SamtoolsStatsFiles
}
