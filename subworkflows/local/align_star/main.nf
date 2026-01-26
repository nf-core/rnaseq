//
// Alignment with STAR and optional UMI deduplication
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN_IGENOMES                       } from '../../../modules/local/star_align_igenomes'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'
include { BAM_DEDUP_UMI                             } from '../../nf-core/bam_dedup_umi'


//
// Function that parses and returns the alignment rate from the STAR log output
//
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
    star_ignore_sjdbgtf     // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    seq_platform            // string : sequencing platform
    seq_center              // string : sequencing center
    is_aws_igenome          // boolean: whether the genome files are from AWS iGenomes
    fasta                   // channel: /path/to/fasta
    use_sentieon_star       // boolean: whether star alignment is accelerated with Sentieon
    with_umi                // boolean: whether UMI processing is enabled
    umi_dedup_tool          // string: 'umicollapse' or 'umitools'
    umitools_dedup_stats    // boolean: whether to generate UMI-tools dedup stats
    bam_csi_index           // boolean: whether to generate CSI index
    skip_markduplicates     // boolean: skip marking duplicates
    transcript_fasta        // channel: [ val(meta), [ fasta ] ] - for UMI dedup
    input_genome_bam        // channel: [ val(meta), path(bam) ] - pre-existing genome BAMs to mix (when not using UMI)
    input_genome_bam_index  // channel: [ val(meta), path(bai) ] - pre-existing genome BAM indices
    input_transcriptome_bam // channel: [ val(meta), path(bam) ] - pre-existing transcriptome BAMs to mix

    main:

    ch_versions = channel.empty()

    //
    // Map reads with STAR
    //
    ch_star_out = null
    if (use_sentieon_star) {

        SENTIEON_STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = SENTIEON_STAR_ALIGN
        // SENTIEON_STAR_ALIGN uses topic-based version reporting

    } else if (is_aws_igenome) {

        STAR_ALIGN_IGENOMES(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = STAR_ALIGN_IGENOMES
        ch_versions = ch_versions.mix(STAR_ALIGN_IGENOMES.out.versions.first())

    } else {

        STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        ch_star_out = STAR_ALIGN
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    }

    ch_orig_bam = ch_star_out.out.bam
    ch_log_final = ch_star_out.out.log_final
    ch_log_out = ch_star_out.out.log_out
    ch_log_progress = ch_star_out.out.log_progress
    ch_bam_sorted = ch_star_out.out.bam_sorted
    ch_bam_transcript = ch_star_out.out.bam_transcript
    ch_fastq = ch_star_out.out.fastq
    ch_tab = ch_star_out.out.tab
    ch_percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    ch_genome_bam       = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_genome_bam_index = bam_csi_index ? BAM_SORT_STATS_SAMTOOLS.out.csi : BAM_SORT_STATS_SAMTOOLS.out.bai
    ch_transcriptome_bam = ch_bam_transcript

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

    // Initialize multiqc_files with STAR log (always added)
    ch_multiqc_files = ch_log_final.collect{ tuple -> tuple[1] }

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //
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

        ch_genome_bam        = BAM_DEDUP_UMI.out.bam
        ch_genome_bam_index  = BAM_DEDUP_UMI.out.bai
        ch_transcriptome_bam = BAM_DEDUP_UMI.out.transcriptome_bam
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
        ch_genome_bam        = input_genome_bam.mix(ch_genome_bam)
        ch_genome_bam_index  = input_genome_bam_index.mix(ch_genome_bam_index)
        ch_transcriptome_bam = input_transcriptome_bam.mix(ch_transcriptome_bam)

        if (skip_markduplicates) {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise. If mark duplicates
            // will run, those stats will be added later instead to avoid
            // duplicate flagstat files in MultiQC.
            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect{ tuple -> tuple[1] })
                .mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect{ tuple -> tuple[1] })
                .mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{ tuple -> tuple[1] })
        }
    }

    emit:
    bam             = ch_genome_bam                          // channel: [ val(meta), path(bam) ]
    bai             = ch_genome_bam_index                    // channel: [ val(meta), path(bai) ]
    csi             = BAM_SORT_STATS_SAMTOOLS.out.csi        // channel: [ val(meta), path(csi) ]
    orig_bam        = ch_orig_bam                            // channel: [ val(meta), path(bam) ] - original aligned BAM before sort/dedup
    bam_transcript  = ch_transcriptome_bam                   // channel: [ val(meta), path(bam) ] - transcriptome BAM (deduplicated if UMI)
    log_final       = ch_log_final                           // channel: [ val(meta), log_final      ]
    log_out         = ch_log_out                             // channel: [ val(meta), log_out        ]
    log_progress    = ch_log_progress                        // channel: [ val(meta), log_progress   ]
    bam_sorted      = ch_bam_sorted                          // channel: [ val(meta), bam_sorted     ]
    fastq           = ch_fastq                               // channel: [ val(meta), fastq          ]
    tab             = ch_tab                                 // channel: [ val(meta), tab            ]
    stats           = BAM_SORT_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat        = BAM_SORT_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats        = BAM_SORT_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    percent_mapped  = ch_percent_mapped                      // channel: [ val(meta), percent_mapped ]

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
