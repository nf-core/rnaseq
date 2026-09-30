// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../../nf-core/bam_stats_samtools/types'

record StarLogs {
    log_final:    Path
    log_out:      Path
    log_progress: Path
    tab:          List<Path>
}

// Emitted by STAR_ALIGN, SENTIEON_STARALIGN and PARABRICKS_RNAFQ2BAM alike.
// orig_bai, qc_metrics and duplicate_metrics are only set by Parabricks.
record StarAlignResult {
    id:                 String
    meta:               Map
    orig_bam:           List<Path>
    bam_sorted:         Path?
    bam_sorted_aligned: Path?
    bam_unsorted:       Path?
    transcriptome_bam:  Path?
    unmapped:           List<Path>
    sam:                Path?
    junction:           Path?
    spl_junc_tab:       Path?
    read_per_gene_tab:  Path?
    wig:                List<Path>
    bedgraph:           List<Path>
    orig_bai:           Path?
    qc_metrics:         Path?
    duplicate_metrics:  Path?
    star:               StarLogs
}

// StarAlignResult plus aligner, percent_mapped and the sorted BAM with its samtools stats.
record StarAligned {
    id:                 String
    meta:               Map
    aligner:            String
    orig_bam:           List<Path>
    bam_sorted:         Path?
    bam_sorted_aligned: Path?
    bam_unsorted:       Path?
    transcriptome_bam:  Path?
    unmapped:           List<Path>
    sam:                Path?
    junction:           Path?
    spl_junc_tab:       Path?
    read_per_gene_tab:  Path?
    wig:                List<Path>
    bedgraph:           List<Path>
    orig_bai:           Path?
    qc_metrics:         Path?
    duplicate_metrics:  Path?
    percent_mapped:     Float
    star:               StarLogs
    bam:                Path
    bai:                Path
    samtools:           SamtoolsStatsFiles
}
