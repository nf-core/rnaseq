//
// Alignment with STAR
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { PARABRICKS_RNAFQ2BAM as PARABRICKS_RNA_FQ2BAM } from '../../../modules/nf-core/parabricks/rnafq2bam/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'
include { StarAligned                               } from './types'


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
    reads                // channel: [ val(meta), [ reads ] ]
    index                // channel: [ val(meta), [ index ] ]
    gtf                  // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf  // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    use_sentieon_star    // boolean: whether star alignment is accelerated with Sentieon
    use_parabricks_star  // boolean: whether star alignment (and mark duplicates) is accelerated with Parabricks
    skip_markduplicates  // boolean: whether to skip marking duplicates

    main:

    //
    // Map reads with STAR
    //
    // Single-end reads arrive as a lone path; the aligner processes take a list.
    ch_reads = reads.map { meta, fastqs -> [ meta, [ fastqs ].flatten() ] }
    ch_star_out = null
    if (use_sentieon_star) {

        SENTIEON_STAR_ALIGN(ch_reads, index, gtf, star_ignore_sjdbgtf)
        ch_star_out = SENTIEON_STAR_ALIGN
        // SENTIEON_STAR_ALIGN uses topic-based version reporting

    } else if (use_parabricks_star) {

        PARABRICKS_RNA_FQ2BAM(ch_reads, fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }, index, true, !skip_markduplicates)
        ch_star_out = PARABRICKS_RNA_FQ2BAM

    } else {

        STAR_ALIGN(ch_reads, index, gtf, star_ignore_sjdbgtf)
        ch_star_out = STAR_ALIGN

    }

    // A run that produced no BAM drops out of the downstream channels.
    ch_star = ch_star_out.out.filter { r -> r.orig_bam }
    ch_orig_bam = ch_star.map { r -> [ r.meta, r.orig_bam ] }
    ch_log_final = ch_star.map { r -> [ r.meta, r.star.log_final ] }
    ch_log_out = ch_star.map { r -> [ r.meta, r.star.log_out ] }
    ch_log_progress = ch_star.map { r -> [ r.meta, r.star.log_progress ] }
    ch_bam_sorted = ch_star.filter { r -> r.bam_sorted }.map { r -> [ r.meta, r.bam_sorted ] }
    ch_bam_transcript = ch_star.filter { r -> r.transcriptome_bam }.map { r -> [ r.meta, r.transcriptome_bam ] }
    ch_fastq = ch_star.filter { r -> r.unmapped }.map { r -> [ r.meta, r.unmapped ] }
    ch_tab = ch_star.filter { r -> r.star.tab }.map { r -> [ r.meta, r.star.tab ] }
    ch_percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta_fai)

    ch_results = ch_star
        .map { r -> r + record(aligner: 'star', percent_mapped: getStarPercentMapped(params, r.star.log_final)) }
        .join(BAM_SORT_STATS_SAMTOOLS.out.results, by: 'id')

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
    index = BAM_SORT_STATS_SAMTOOLS.out.index       // channel: [ val(meta), [ index ] ]
    stats = BAM_SORT_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped              // channel: [ val(meta), percent_mapped ]
    results = ch_results                            // channel: StarAligned
}
