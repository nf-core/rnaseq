// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'

// umi_tools writes all three stats tables together, so the tsv fields are either all set or all null.
record UmitoolsDedupBam {
    id:                   String
    meta:                 Map
    bam:                  Path
    bai:                  Path
    log:                  Path
    tsv_edit_distance:    Path?
    tsv_per_umi:          Path?
    tsv_umi_per_position: Path?
    samtools:             SamtoolsStatsFiles
}

record UmitoolsDedupStats {
    edit_distance:    Path
    per_umi:          Path
    umi_per_position: Path
}
