// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'
include { UmitoolsDedupStats } from '../bam_dedup_stats_samtools_umitools/types'

record UmiDedupTranscriptome {
    transcriptome_bam:      Path
    dedup_bam:              Path
    sorted_bam:             Path
    sorted_bam_index:       Path
    filtered_bam:           Path?
    samtools:               SamtoolsStatsFiles
    tsv:                    UmitoolsDedupStats?
    // Coordinate-sorted, pre-dedup form (built once, before the umitools/umicollapse branch).
    coord_sorted_bam:       Path?
    coord_sorted_bam_index: Path?
    coord_sorted_samtools:  SamtoolsStatsFiles?
}

record UmiDedupBam {
    id:                       String
    meta:                     Map
    bam:                      Path
    bai:                      Path
    samtools:                 SamtoolsStatsFiles
    genomic_dedup_log:        Path
    transcriptomic_dedup_log: Path?
    prepare_for_rsem_log:     Path?
    transcriptome_bam:        Path?
    transcriptome:            UmiDedupTranscriptome?
    tsv:                      UmitoolsDedupStats?
}
