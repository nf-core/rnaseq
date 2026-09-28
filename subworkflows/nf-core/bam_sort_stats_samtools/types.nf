// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

record SortedBam {
    id:       String
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}
