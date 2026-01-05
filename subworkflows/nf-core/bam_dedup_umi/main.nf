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
    ch_versions = Channel.empty()

    if (umi_dedup_tool != "umicollapse" && umi_dedup_tool != "umitools"){
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME (
            ch_genome_bam
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME
        ch_dedup_log = UMI_DEDUP_GENOME.out.dedup_stats

    } else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
            ch_genome_bam,
            umitools_dedup_stats
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME
        ch_dedup_log = UMI_DEDUP_GENOME.out.deduplog
    }

    // Co-ordinate sort, index and run stats on transcriptome BAM. This takes
    // some preparation- we have to coordinate sort the BAM, run the
    // deduplication, then restore name sorting and run a script from umitools
    // to prepare for rsem or salmon

    // 1. Coordinate sort

    BAM_SORT_STATS_SAMTOOLS (
        ch_transcriptome_bam,
        ch_transcript_fasta
    )
    ch_sorted_transcriptome_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        .join(BAM_SORT_STATS_SAMTOOLS.out.bai)

    // 2. Transcriptome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME (
            ch_sorted_transcriptome_bam
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME
        ch_dedup_log = ch_dedup_log.mix(UMI_DEDUP_TRANSCRIPTOME.out.dedup_stats)

    } else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
            ch_sorted_transcriptome_bam,
            umitools_dedup_stats
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME
        ch_dedup_log = ch_dedup_log.mix(UMI_DEDUP_TRANSCRIPTOME.out.deduplog)
    }

    // 3. Restore name sorting
    SAMTOOLS_SORT (
        UMI_DEDUP_TRANSCRIPTOME.out.bam,
        ch_fasta
    )

    // 4. Run prepare_for_rsem.py on paired-end BAM files
    // This fixes paired-end reads in name sorted BAM files
    // See: https://github.com/nf-core/rnaseq/issues/828
    ended_transcriptome_dedup_bam = SAMTOOLS_SORT.out.bam
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

    ch_dedup_transcriptome_bam = ended_transcriptome_dedup_bam.single_end
        .mix(UMITOOLS_PREPAREFORRSEM.out.bam)

    // Collect files useful for MultiQC into one helpful emission. Don't
    // automatically add transcriptome stats- difficult to separate in multiqc
    // without a bit more work

    ch_multiqc_files = ch_dedup_log
        .mix(UMI_DEDUP_GENOME.out.stats)
        .mix(UMI_DEDUP_GENOME.out.flagstat)
        .mix(UMI_DEDUP_GENOME.out.idxstats)
        .transpose()
        .map{it[1]}

    // Record versions

    ch_versions = UMI_DEDUP_GENOME.out.versions
        .mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
        .mix(UMITOOLS_PREPAREFORRSEM.out.versions)

    emit:
    bam                = UMI_DEDUP_GENOME.out.bam                                                // channel: [ val(meta), path(bam) ]
    bai                = bam_csi_index ? UMI_DEDUP_GENOME.out.csi : UMI_DEDUP_GENOME.out.bai     // channel: [ val(meta), path(bai) ]
    dedup_log          = ch_dedup_log                                                            // channel: [ val(meta), path(log) ]
    stats              = UMI_DEDUP_GENOME.out.stats.mix(UMI_DEDUP_TRANSCRIPTOME.out.stats)       // channel: [ val(meta), path(stats)]
    flagstat           = UMI_DEDUP_GENOME.out.flagstat.mix(UMI_DEDUP_TRANSCRIPTOME.out.flagstat) // channel: [ val(meta), path(flagstat)]
    idxstats           = UMI_DEDUP_GENOME.out.idxstats.mix(UMI_DEDUP_TRANSCRIPTOME.out.idxstats) // channel: [ val(meta), path(idxstats)]
    multiqc_files      = ch_multiqc_files                                                        // channel: file
    transcriptome_bam  = ch_dedup_transcriptome_bam                                              // channel: [ val(meta), path(bam) ]
    versions            = ch_versions                                                            // channel: [ path(versions.yml) ]
}
