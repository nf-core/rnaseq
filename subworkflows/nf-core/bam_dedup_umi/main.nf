nextflow.preview.types = true

//
// BAM deduplication with UMI processing
//

include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME       } from '../bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME        } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME              } from '../bam_dedup_stats_samtools_umitools'
include { BAM_SORT_STATS_SAMTOOLS                                                                    } from '../bam_sort_stats_samtools'

include { UMITOOLS_PREPAREFORRSEM  } from '../../../modules/nf-core/umitools/prepareforrsem'
include { SAMTOOLS_SORT            } from '../../../modules/nf-core/samtools/sort/main'

record UmiDedupResult {
    meta:                        Map
    bam:                         Path
    bai:                         Path
    genomic_dedup_log:           Path?
    transcriptomic_dedup_log:    Path?
    prepare_for_rsem_log:        Path?
    transcriptome_bam:           Path?
    transcriptome_dedup_bam:     Path?
    transcriptome_sorted_bam:    Path?
    transcriptome_sorted_bam_bai: Path?
    transcriptome_filtered_bam:  Path?
    genome_stats:                Path?
    genome_flagstat:             Path?
    genome_idxstats:             Path?
    transcriptome_stats:         Path?
    transcriptome_flagstat:      Path?
    transcriptome_idxstats:      Path?
    tsv_edit_distance:           Path?
    tsv_per_umi:                 Path?
    tsv_umi_per_position:        Path?
}

workflow BAM_DEDUP_UMI {
    take:
    ch_genome_bam         // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta              // channel: [ val(meta), path(fasta) ]
    umi_dedup_tool        // string: 'umicollapse' or 'umitools'
    umitools_dedup_stats  // boolean: whether to generate UMI-tools dedup stats
    bam_csi_index         // boolean: whether to generate CSI index
    ch_transcriptome_bam  // channel: [ val(meta), path(bam) ]
    ch_transcript_fasta   // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = channel.empty()
    ch_tsv_edit_distance    = channel.empty()
    ch_tsv_per_umi          = channel.empty()
    ch_tsv_umi_per_position = channel.empty()
    ch_genomic_dedup_log       = channel.empty()
    ch_transcriptomic_dedup_log = channel.empty()

    if (umi_dedup_tool != "umicollapse" && umi_dedup_tool != "umitools"){
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME (
            ch_genome_bam
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME
        ch_genomic_dedup_log = UMI_DEDUP_GENOME.out.dedup_stats

    } else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
            ch_genome_bam,
            umitools_dedup_stats
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME
        ch_genomic_dedup_log = UMI_DEDUP_GENOME.out.deduplog
        ch_tsv_edit_distance    = UMI_DEDUP_GENOME.out.tsv_edit_distance
        ch_tsv_per_umi          = UMI_DEDUP_GENOME.out.tsv_per_umi
        ch_tsv_umi_per_position = UMI_DEDUP_GENOME.out.tsv_umi_per_position
    }

    BAM_SORT_STATS_SAMTOOLS (
        ch_transcriptome_bam,
        ch_transcript_fasta
    )

    // Record field access: extract [meta, bam, bai] from SamtoolsResult record
    ch_sorted_transcriptome_bam = BAM_SORT_STATS_SAMTOOLS.out.result
        .map { r -> [r.meta, r.bam, r.bai] }

    // Transcriptome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME (
            ch_sorted_transcriptome_bam
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME
        ch_transcriptomic_dedup_log = UMI_DEDUP_TRANSCRIPTOME.out.dedup_stats

    } else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
            ch_sorted_transcriptome_bam,
            umitools_dedup_stats
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME
        ch_transcriptomic_dedup_log = UMI_DEDUP_TRANSCRIPTOME.out.deduplog
        ch_tsv_edit_distance    = ch_tsv_edit_distance.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_edit_distance)
        ch_tsv_per_umi          = ch_tsv_per_umi.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_per_umi)
        ch_tsv_umi_per_position = ch_tsv_umi_per_position.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_umi_per_position)
    }

    // Restore name sorting
    SAMTOOLS_SORT (
        UMI_DEDUP_TRANSCRIPTOME.out.bam,
        ch_fasta,
        ''
    )

    // Extract bam from SamtoolsSortResult record
    ch_namesorted_bam = SAMTOOLS_SORT.out.map { r -> [r.meta, r.bam] }

    // Fix paired-end reads in name sorted BAM files
    ended_transcriptome_dedup_bam = ch_namesorted_bam
        .branch {
            meta, bam ->
                single_end: meta.single_end
                    return [ meta, bam ]
                paired_end: !meta.single_end
                    return [ meta, bam ]
        }

    UMITOOLS_PREPAREFORRSEM (
        ended_transcriptome_dedup_bam.paired_end
            .map { meta, bam -> [ meta, bam, [] ] }
    )

    // Extract fields from PrepareForRsemResult record
    ch_rsem_bam = UMITOOLS_PREPAREFORRSEM.out.map { r -> [r.meta, r.bam] }
    ch_rsem_log = UMITOOLS_PREPAREFORRSEM.out.map { r -> [r.meta, r.log] }

    ch_dedup_transcriptome_bam = ended_transcriptome_dedup_bam.single_end
        .mix(ch_rsem_bam)

    ch_multiqc_files = ch_genomic_dedup_log
        .mix(UMI_DEDUP_GENOME.out.stats)
        .mix(UMI_DEDUP_GENOME.out.flagstat)
        .mix(UMI_DEDUP_GENOME.out.idxstats)
        .transpose()
        .map{ item -> item[1] }

    ch_versions = UMI_DEDUP_GENOME.out.versions
        .mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // Join all per-sample outputs by meta to construct the record.
    // Genome and transcriptome stats are kept as separate fields
    // (rather than mixed) to enable clean join-by-meta construction.
    ch_bai = bam_csi_index ? UMI_DEDUP_GENOME.out.csi : UMI_DEDUP_GENOME.out.bai

    emit:
    result = UMI_DEDUP_GENOME.out.bam
        .join(ch_bai, by: [0])
        .join(ch_genomic_dedup_log,        by: [0], remainder: true)
        .join(ch_transcriptomic_dedup_log, by: [0], remainder: true)
        .join(ch_rsem_log, by: [0], remainder: true)
        .join(ch_dedup_transcriptome_bam,  by: [0], remainder: true)
        .join(UMI_DEDUP_TRANSCRIPTOME.out.bam, by: [0], remainder: true)
        .join(ch_namesorted_bam,       by: [0], remainder: true)
        .join(UMI_DEDUP_TRANSCRIPTOME.out.bai, by: [0], remainder: true)
        .join(ch_rsem_bam, by: [0], remainder: true)
        .join(UMI_DEDUP_GENOME.out.stats,    by: [0], remainder: true)
        .join(UMI_DEDUP_GENOME.out.flagstat,  by: [0], remainder: true)
        .join(UMI_DEDUP_GENOME.out.idxstats,  by: [0], remainder: true)
        .join(UMI_DEDUP_TRANSCRIPTOME.out.stats,    by: [0], remainder: true)
        .join(UMI_DEDUP_TRANSCRIPTOME.out.flagstat,  by: [0], remainder: true)
        .join(UMI_DEDUP_TRANSCRIPTOME.out.idxstats,  by: [0], remainder: true)
        .join(ch_tsv_edit_distance,    by: [0], remainder: true)
        .join(ch_tsv_per_umi,          by: [0], remainder: true)
        .join(ch_tsv_umi_per_position, by: [0], remainder: true)
        .map { meta, bam, bai, g_log, t_log, rsem_log, t_bam, t_dedup_bam, t_sorted_bam, t_sorted_bai, t_filtered_bam,
               g_stats, g_flagstat, g_idxstats, t_stats, t_flagstat, t_idxstats,
               tsv_ed, tsv_pu, tsv_upp ->
            record(
                meta: meta,
                bam: bam, bai: bai,
                genomic_dedup_log: g_log, transcriptomic_dedup_log: t_log,
                prepare_for_rsem_log: rsem_log,
                transcriptome_bam: t_bam, transcriptome_dedup_bam: t_dedup_bam,
                transcriptome_sorted_bam: t_sorted_bam, transcriptome_sorted_bam_bai: t_sorted_bai,
                transcriptome_filtered_bam: t_filtered_bam,
                genome_stats: g_stats, genome_flagstat: g_flagstat, genome_idxstats: g_idxstats,
                transcriptome_stats: t_stats, transcriptome_flagstat: t_flagstat, transcriptome_idxstats: t_idxstats,
                tsv_edit_distance: tsv_ed, tsv_per_umi: tsv_pu, tsv_umi_per_position: tsv_upp
            )
        }
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
