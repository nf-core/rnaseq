nextflow.enable.types = true

//
// BAM deduplication with UMI processing
//

include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME       } from '../bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME        } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME              } from '../bam_dedup_stats_samtools_umitools'
include { BAM_SORT_STATS_SAMTOOLS                                                                    } from '../bam_sort_stats_samtools'

include { UMITOOLS_PREPAREFORRSEM                                                                    } from '../../../modules/nf-core/umitools/prepareforrsem'
include { SAMTOOLS_SORT                                                                              } from '../../../modules/nf-core/samtools/sort/main'
include { Bam; UmiDedupBam } from '../../../modules/nf-core/types'
include { SamtoolsSortResult } from '../../../modules/nf-core/samtools/sort/main'
include { UmitoolsPrepareforrsemResult } from '../../../modules/nf-core/umitools/prepareforrsem/main'

workflow BAM_DEDUP_UMI {
    take:
    ch_genome_bam: Channel<Bam>
    fasta: Value<Path?>
    fai: Value<Path?>
    umi_dedup_tool: String // 'umicollapse' or 'umitools'
    umitools_dedup_stats: Boolean // whether to generate UMI-tools dedup stats
    ch_transcriptome_bam: Channel<Bam> // records with transcriptome_bam set
    transcript_fasta: Value<Path?>
    umitools_dedup_primary_only: Boolean // whether to filter to primary alignments before dedup
    tool_args: Record // the samtools index options, forwarded
    umi_grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:
    if (umi_dedup_tool != "umicollapse" && umi_dedup_tool != "umitools") {
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    // Both tools are reduced to the same record shape; only umitools produces stats tables (tsv).

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        ch_genome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME(ch_genome_bam, tool_args, umi_grouping_method, umi_separator)
            .map { r -> r + record(tsv: null) }
    }
    else {
        ch_genome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME(
            ch_genome_bam,
            umitools_dedup_stats,
            umitools_dedup_primary_only,
            tool_args,
            umi_grouping_method,
            umi_separator,
        )
            .map { r ->
                record(
                    id: r.id, meta: r.meta, bam: r.bam, bai: r.bai, log: r.log, samtools: r.samtools,
                    tsv: r.tsv_edit_distance != null
                        ? record(edit_distance: r.tsv_edit_distance, per_umi: r.tsv_per_umi, umi_per_position: r.tsv_umi_per_position)
                        : null
                )
            }
    }

    // Co-ordinate sort, index and run stats on transcriptome BAM. This takes
    // some preparation- we have to coordinate sort the BAM, run the
    // deduplication, then restore name sorting and run a script from umitools
    // to prepare for rsem or salmon

    // 1. Coordinate sort
    def ch_coord_sorted: Channel<Bam> = BAM_SORT_STATS_SAMTOOLS(
        ch_transcriptome_bam.map { r -> record(id: r.id, meta: r.meta, raw_bams: [r.transcriptome_bam]) },
        transcript_fasta,
        null,
        tool_args
    )

    // 2. Transcriptome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        ch_transcriptome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME(ch_coord_sorted, tool_args, umi_grouping_method, umi_separator)
            .map { r -> r + record(tsv: null) }
    }
    else {
        ch_transcriptome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME(
            ch_coord_sorted,
            umitools_dedup_stats,
            umitools_dedup_primary_only,
            tool_args,
            umi_grouping_method,
            umi_separator,
        )
            .map { r ->
                record(
                    id: r.id, meta: r.meta, bam: r.bam, bai: r.bai, log: r.log, samtools: r.samtools,
                    tsv: r.tsv_edit_distance != null
                        ? record(edit_distance: r.tsv_edit_distance, per_umi: r.tsv_per_umi, umi_per_position: r.tsv_umi_per_position)
                        : null
                )
            }
    }

    // 3. Restore name sorting
    def ch_name_sorted: Channel<SamtoolsSortResult> = SAMTOOLS_SORT(
        ch_transcriptome_dedup.map { r -> record(id: r.id, meta: r.meta, raw_bams: [r.bam]) },
        fasta,
        fai,
        '',
    )
        .filter { r -> r.bam != null }

    // 4. Run prepare_for_rsem.py on paired-end BAM files
    // This fixes paired-end reads in name sorted BAM files
    // See: https://github.com/nf-core/rnaseq/issues/828
    def ch_prepared: Channel<UmitoolsPrepareforrsemResult> = UMITOOLS_PREPAREFORRSEM(ch_name_sorted.filter { r -> !r.meta.single_end })

    // Only paired-end samples pass through UMITOOLS_PREPAREFORRSEM, so the remainder
    // join leaves `prepared` null for single-end samples.
    ch_transcriptome_part = ch_transcriptome_bam
        .join(ch_coord_sorted.map { r -> record(id: r.id, coord_sorted: r) }, by: 'id')
        .join(ch_transcriptome_dedup.map { r -> record(id: r.id, dedup: r) }, by: 'id')
        .join(ch_name_sorted.map { r -> record(id: r.id, name_sorted: r) }, by: 'id')
        .join(ch_prepared.map { r -> record(id: r.id, prepared: r) }, by: 'id', remainder: true)
        .map { r ->
            record(
                id:                       r.id,
                transcriptomic_dedup_log: r.dedup.log,
                prepare_for_rsem_log:     r.prepared?.log,
                transcriptome_bam:        r.prepared?.bam ?: r.name_sorted.bam,
                transcriptome:            record(
                    transcriptome_bam:      r.transcriptome_bam,
                    dedup_bam:              r.dedup.bam,
                    sorted_bam:             r.name_sorted.bam,
                    sorted_bam_index:       r.dedup.bai,
                    filtered_bam:           r.prepared?.bam,
                    samtools:               r.dedup.samtools,
                    tsv:                    r.dedup.tsv,
                    coord_sorted_bam:       r.coord_sorted.bam,
                    coord_sorted_bam_index: r.coord_sorted.bai,
                    coord_sorted_samtools:  r.coord_sorted.samtools
                )
            )
        }

    // The transcriptome side is absent when there is no transcriptome BAM (e.g. HISAT2).
    ch_results = ch_genome_dedup
        .map { r ->
            record(
                id:                r.id,
                meta:              r.meta,
                bam:               r.bam,
                bai:               r.bai,
                samtools:          r.samtools,
                genomic_dedup_log: r.log,
                tsv:               r.tsv
            )
        }
        .join(ch_transcriptome_part, by: 'id', remainder: true)

    emit:
    ch_results
}
