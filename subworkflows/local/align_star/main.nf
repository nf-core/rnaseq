//
// Alignment with STAR
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN_IGENOMES                       } from '../../../modules/local/star_align_igenomes'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'


//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(params, align_log) {
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
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    fasta               // channel: /path/to/fasta
    use_sentieon_star   // boolean: whether star alignment is accelerated with Sentieon

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
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

    ch_orig_bam = ch_star_out.out.bam
    ch_log_final = ch_star_out.out.log_final
    ch_log_out = ch_star_out.out.log_out
    ch_log_progress = ch_star_out.out.log_progress
    ch_bam_sorted = ch_star_out.out.bam_sorted
    ch_bam_transcript = ch_star_out.out.bam_transcript
    ch_fastq = ch_star_out.out.fastq
    ch_tab = ch_star_out.out.tab
    ch_versions = ch_versions.mix(ch_star_out.out.versions.first())
    ch_percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    orig_bam = ch_orig_bam                          // channel: [ val(meta), bam            ]
    log_final = ch_log_final                        // channel: [ val(meta), log_final      ]
    log_out = ch_log_out                            // channel: [ val(meta), log_out        ]
    log_progress = ch_log_progress                  // channel: [ val(meta), log_progress   ]
    bam_sorted = ch_bam_sorted                      // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq = ch_fastq                                // channel: [ val(meta), fastq          ]
    tab = ch_tab                                    // channel: [ val(meta), tab            ]
    bam = BAM_SORT_STATS_SAMTOOLS.out.bam           // channel: [ val(meta), [ bam ] ]
    bai = BAM_SORT_STATS_SAMTOOLS.out.bai           // channel: [ val(meta), [ bai ] ]
    csi = BAM_SORT_STATS_SAMTOOLS.out.csi           // channel: [ val(meta), [ csi ] ]
    stats = BAM_SORT_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped              // channel: [ val(meta), percent_mapped ]
    versions = ch_versions                          // channel: [ versions.yml ]
}
