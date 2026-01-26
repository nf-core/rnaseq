//
// Alignment with HISAT2 and optional UMI deduplication
//

include { FASTQ_ALIGN_HISAT2 } from '../../nf-core/fastq_align_hisat2'
include { BAM_DEDUP_UMI      } from '../../nf-core/bam_dedup_umi'

workflow ALIGN_HISAT2 {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    index                   // channel: [ val(meta), [ index ] ]
    splicesites             // channel: [ val(meta), [ splicesites ] ]
    fasta                   // channel: [ val(meta), [ fasta ] ]
    with_umi                // boolean: whether UMI processing is enabled
    umi_dedup_tool          // string: 'umicollapse' or 'umitools'
    umitools_dedup_stats    // boolean: whether to generate UMI-tools dedup stats
    bam_csi_index           // boolean: whether to generate CSI index
    skip_markduplicates     // boolean: skip marking duplicates
    transcriptome_bam       // channel: [ val(meta), path(bam) ] - for UMI dedup
    transcript_fasta        // channel: [ val(meta), [ fasta ] ] - for UMI dedup
    input_genome_bam        // channel: [ val(meta), path(bam) ] - pre-existing genome BAMs to mix (when not using UMI)
    input_genome_bam_index  // channel: [ val(meta), path(bai) ] - pre-existing genome BAM indices

    main:
    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Align reads with HISAT2
    //
    FASTQ_ALIGN_HISAT2 (
        reads,
        index,
        splicesites,
        fasta
    )

    ch_genome_bam       = FASTQ_ALIGN_HISAT2.out.bam
    ch_genome_bam_index = bam_csi_index ? FASTQ_ALIGN_HISAT2.out.csi : FASTQ_ALIGN_HISAT2.out.bai
    ch_versions         = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

    // Initialize UMI output channels as empty - will be populated if with_umi is true
    ch_umi_genomic_dedup_log        = channel.empty()
    ch_umi_transcriptomic_dedup_log = channel.empty()
    ch_umi_prepare_for_rsem_log     = channel.empty()
    ch_umi_transcriptome_dedup_bam      = channel.empty()
    ch_umi_transcriptome_sorted_bam     = channel.empty()
    ch_umi_transcriptome_sorted_bam_bai = channel.empty()
    ch_umi_transcriptome_filtered_bam   = channel.empty()
    ch_umi_dedup_stats    = channel.empty()
    ch_umi_dedup_bam      = channel.empty()
    ch_umi_dedup_bai      = channel.empty()
    ch_umi_dedup_flagstat = channel.empty()
    ch_umi_dedup_idxstats = channel.empty()
    ch_umi_dedup_tsv_edit_distance    = channel.empty()
    ch_umi_dedup_tsv_per_umi          = channel.empty()
    ch_umi_dedup_tsv_umi_per_position = channel.empty()

    // Initialize multiqc_files with HISAT2 summary (always added)
    ch_multiqc_files = FASTQ_ALIGN_HISAT2.out.summary.collect{ tuple -> tuple[1] }

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //
    if (with_umi) {
        def ch_bam_for_dedup = input_genome_bam.mix(ch_genome_bam)
        def ch_bai_for_dedup = input_genome_bam_index.mix(ch_genome_bam_index)

        BAM_DEDUP_UMI(
            ch_bam_for_dedup.join(ch_bai_for_dedup, by: [0]),
            fasta,
            umi_dedup_tool,
            umitools_dedup_stats,
            bam_csi_index,
            transcriptome_bam,
            transcript_fasta
        )

        ch_genome_bam       = BAM_DEDUP_UMI.out.bam
        ch_genome_bam_index = BAM_DEDUP_UMI.out.bai
        ch_umi_genomic_dedup_log        = BAM_DEDUP_UMI.out.genomic_dedup_log
        ch_umi_transcriptomic_dedup_log = BAM_DEDUP_UMI.out.transcriptomic_dedup_log
        ch_umi_prepare_for_rsem_log     = BAM_DEDUP_UMI.out.prepare_for_rsem_log
        ch_umi_transcriptome_dedup_bam      = BAM_DEDUP_UMI.out.transcriptome_dedup_bam
        ch_umi_transcriptome_sorted_bam     = BAM_DEDUP_UMI.out.transcriptome_sorted_bam
        ch_umi_transcriptome_sorted_bam_bai = BAM_DEDUP_UMI.out.transcriptome_sorted_bam_bai
        ch_umi_transcriptome_filtered_bam   = BAM_DEDUP_UMI.out.transcriptome_filtered_bam
        ch_umi_dedup_stats    = BAM_DEDUP_UMI.out.stats
        ch_umi_dedup_bam      = BAM_DEDUP_UMI.out.bam
        ch_umi_dedup_bai      = BAM_DEDUP_UMI.out.bai
        ch_umi_dedup_flagstat = BAM_DEDUP_UMI.out.flagstat
        ch_umi_dedup_idxstats = BAM_DEDUP_UMI.out.idxstats
        ch_umi_dedup_tsv_edit_distance    = BAM_DEDUP_UMI.out.tsv_edit_distance
        ch_umi_dedup_tsv_per_umi          = BAM_DEDUP_UMI.out.tsv_per_umi
        ch_umi_dedup_tsv_umi_per_position = BAM_DEDUP_UMI.out.tsv_umi_per_position
        ch_versions      = ch_versions.mix(BAM_DEDUP_UMI.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_UMI.out.multiqc_files)

    } else {
        ch_genome_bam       = input_genome_bam.mix(ch_genome_bam)
        ch_genome_bam_index = input_genome_bam_index.mix(ch_genome_bam_index)

        if (skip_markduplicates) {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise. If mark duplicates
            // will run, those stats will be added later instead to avoid
            // duplicate flagstat files in MultiQC.
            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{ tuple -> tuple[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{ tuple -> tuple[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{ tuple -> tuple[1] })
        }
    }

    emit:
    bam             = ch_genome_bam                          // channel: [ val(meta), path(bam) ]
    bai             = ch_genome_bam_index                    // channel: [ val(meta), path(bai) ]
    orig_bam        = FASTQ_ALIGN_HISAT2.out.bam             // channel: [ val(meta), path(bam) ] - original aligned BAM before dedup
    unaligned       = FASTQ_ALIGN_HISAT2.out.fastq           // channel: [ val(meta), path(fastq) ]
    summary         = FASTQ_ALIGN_HISAT2.out.summary         // channel: [ val(meta), path(summary) ]
    stats           = FASTQ_ALIGN_HISAT2.out.stats           // channel: [ val(meta), path(stats) ]
    flagstat        = FASTQ_ALIGN_HISAT2.out.flagstat        // channel: [ val(meta), path(flagstat) ]
    idxstats        = FASTQ_ALIGN_HISAT2.out.idxstats        // channel: [ val(meta), path(idxstats) ]

    // UMI dedup outputs
    umi_genomic_dedup_log        = ch_umi_genomic_dedup_log
    umi_transcriptomic_dedup_log = ch_umi_transcriptomic_dedup_log
    umi_prepare_for_rsem_log     = ch_umi_prepare_for_rsem_log
    umi_transcriptome_dedup_bam      = ch_umi_transcriptome_dedup_bam
    umi_transcriptome_sorted_bam     = ch_umi_transcriptome_sorted_bam
    umi_transcriptome_sorted_bam_bai = ch_umi_transcriptome_sorted_bam_bai
    umi_transcriptome_filtered_bam   = ch_umi_transcriptome_filtered_bam
    umi_dedup_stats    = ch_umi_dedup_stats
    umi_dedup_bam      = ch_umi_dedup_bam
    umi_dedup_bai      = ch_umi_dedup_bai
    umi_dedup_flagstat = ch_umi_dedup_flagstat
    umi_dedup_idxstats = ch_umi_dedup_idxstats
    umi_dedup_tsv_edit_distance    = ch_umi_dedup_tsv_edit_distance
    umi_dedup_tsv_per_umi          = ch_umi_dedup_tsv_per_umi
    umi_dedup_tsv_umi_per_position = ch_umi_dedup_tsv_umi_per_position

    multiqc_files   = ch_multiqc_files                       // channel: [ path(files) ]
    versions        = ch_versions                            // channel: [ versions.yml ]
}
