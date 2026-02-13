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
    input_genome_bam        // channel: [ val(meta), path(bam) ]
    input_genome_bam_index  // channel: [ val(meta), path(bai) ]

    main:
    ch_versions = channel.empty()

    FASTQ_ALIGN_HISAT2 (
        reads,
        index,
        splicesites,
        fasta
    )

    ch_genome_bam       = FASTQ_ALIGN_HISAT2.out.bam
    ch_genome_bam_index = bam_csi_index ? FASTQ_ALIGN_HISAT2.out.csi : FASTQ_ALIGN_HISAT2.out.bai
    ch_versions         = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

    ch_umi_result = channel.empty()
    ch_multiqc_files = FASTQ_ALIGN_HISAT2.out.summary.collect{ tuple -> tuple[1] }

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

        ch_umi_result       = BAM_DEDUP_UMI.out.result
        ch_genome_bam       = ch_umi_result.map { r -> [r.meta, r.bam] }
        ch_genome_bam_index = ch_umi_result.map { r -> [r.meta, r.bai] }
        ch_versions      = ch_versions.mix(BAM_DEDUP_UMI.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_UMI.out.multiqc_files)

    } else {
        ch_genome_bam       = input_genome_bam.mix(ch_genome_bam)
        ch_genome_bam_index = input_genome_bam_index.mix(ch_genome_bam_index)

        if (skip_markduplicates) {
            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{ tuple -> tuple[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{ tuple -> tuple[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{ tuple -> tuple[1] })
        }
    }

    emit:
    bam             = ch_genome_bam
    bai             = ch_genome_bam_index
    orig_bam        = FASTQ_ALIGN_HISAT2.out.bam
    unaligned       = FASTQ_ALIGN_HISAT2.out.fastq
    summary         = FASTQ_ALIGN_HISAT2.out.summary
    stats           = FASTQ_ALIGN_HISAT2.out.stats
    flagstat        = FASTQ_ALIGN_HISAT2.out.flagstat
    idxstats        = FASTQ_ALIGN_HISAT2.out.idxstats
    umi             = ch_umi_result       // channel: UmiDedupResult records or empty
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
