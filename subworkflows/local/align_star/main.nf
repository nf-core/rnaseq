//
// Alignment with STAR and optional UMI deduplication
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN_IGENOMES                       } from '../../../modules/local/star_align_igenomes'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'
include { BAM_DEDUP_UMI                             } from '../../nf-core/bam_dedup_umi'


def getStarPercentMapped(_params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    return percent_aligned
}

workflow ALIGN_STAR {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    index                   // channel: [ val(meta), [ index ] ]
    gtf                     // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf     // boolean
    seq_platform            // string
    seq_center              // string
    is_aws_igenome          // boolean
    fasta                   // channel: /path/to/fasta
    use_sentieon_star       // boolean
    with_umi                // boolean
    umi_dedup_tool          // string
    umitools_dedup_stats    // boolean
    bam_csi_index           // boolean
    skip_markduplicates     // boolean
    transcript_fasta        // channel: [ val(meta), [ fasta ] ]
    input_genome_bam        // channel: [ val(meta), path(bam) ]
    input_genome_bam_index  // channel: [ val(meta), path(bai) ]
    input_transcriptome_bam // channel: [ val(meta), path(bam) ]

    main:

    ch_versions = channel.empty()

    //
    // Map reads with STAR (all variants emit StarAlignResult record)
    //
    ch_star_out = null
    if (use_sentieon_star) {
        SENTIEON_STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = SENTIEON_STAR_ALIGN
    } else if (is_aws_igenome) {
        STAR_ALIGN_IGENOMES(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = STAR_ALIGN_IGENOMES
    } else {
        STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = STAR_ALIGN
    }

    // Access record fields: .out is channel of StarAlignResult records
    ch_star_result = ch_star_out.out
    ch_orig_bam    = ch_star_result.map { r -> [r.meta, r.bam] }
    ch_percent_mapped = ch_star_result.map { r ->
        [r.meta, getStarPercentMapped(params, r.log_final)]
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // SamtoolsResult record contains bam, bai, csi, stats, flagstat, idxstats
    ch_samtools = BAM_SORT_STATS_SAMTOOLS.out.result

    ch_genome_bam       = ch_samtools.map { r -> [r.meta, r.bam] }
    ch_genome_bam_index = ch_samtools.map { r -> [r.meta, bam_csi_index ? r.csi : r.bai] }
    ch_transcriptome_bam = ch_star_result.map { r -> [r.meta, r.bam_transcript] }

    // UMI dedup record (empty channel when UMI not enabled)
    ch_umi_result = channel.empty()
    ch_multiqc_files = ch_star_result.map { r -> r.log_final }

    if (with_umi) {
        def ch_bam_for_dedup = input_genome_bam.mix(ch_genome_bam)
        def ch_bai_for_dedup = input_genome_bam_index.mix(ch_genome_bam_index)
        def ch_transcriptome_for_dedup = input_transcriptome_bam.mix(ch_transcriptome_bam)

        BAM_DEDUP_UMI(
            ch_bam_for_dedup.join(ch_bai_for_dedup, by: [0]),
            fasta,
            umi_dedup_tool,
            umitools_dedup_stats,
            bam_csi_index,
            ch_transcriptome_for_dedup,
            transcript_fasta
        )

        // UmiDedupResult record replaces 15 individual channel captures
        ch_umi_result        = BAM_DEDUP_UMI.out.result
        ch_genome_bam        = ch_umi_result.map { r -> [r.meta, r.bam] }
        ch_genome_bam_index  = ch_umi_result.map { r -> [r.meta, r.bai] }
        ch_transcriptome_bam = ch_umi_result.map { r -> [r.meta, r.transcriptome_bam] }
        ch_versions      = ch_versions.mix(BAM_DEDUP_UMI.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_UMI.out.multiqc_files)

    } else {
        ch_genome_bam        = input_genome_bam.mix(ch_genome_bam)
        ch_genome_bam_index  = input_genome_bam_index.mix(ch_genome_bam_index)
        ch_transcriptome_bam = input_transcriptome_bam.mix(ch_transcriptome_bam)

        if (skip_markduplicates) {
            ch_multiqc_files = ch_multiqc_files
                .mix(ch_samtools.map { r -> r.stats })
                .mix(ch_samtools.map { r -> r.flagstat })
                .mix(ch_samtools.map { r -> r.idxstats })
        }
    }

    // 35 emits -> 8
    emit:
    bam             = ch_genome_bam
    bai             = ch_genome_bam_index
    bam_transcript  = ch_transcriptome_bam
    star            = ch_star_result      // channel: StarAlignResult records
    samtools        = ch_samtools         // channel: SamtoolsResult records
    umi             = ch_umi_result       // channel: UmiDedupResult records or empty
    percent_mapped  = ch_percent_mapped
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
