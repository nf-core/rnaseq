nextflow.enable.types = true

//
// Alignment with STAR
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { PARABRICKS_RNAFQ2BAM as PARABRICKS_RNA_FQ2BAM } from '../../../modules/nf-core/parabricks/rnafq2bam/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'
include { ReadsInput; StarAlignResult; StarAligned; Bam } from '../../../modules/nf-core/types'


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
    ch_samples: Channel<ReadsInput>
    index: Value<Path>
    gtf: Value<Path?>
    star_ignore_sjdbgtf: Boolean // when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    fasta: Value<Path?>
    fai: Value<Path?>
    use_sentieon_star: Boolean // whether star alignment is accelerated with Sentieon
    use_parabricks_star: Boolean // whether star alignment (and mark duplicates) is accelerated with Parabricks
    skip_markduplicates: Boolean // whether to skip marking duplicates

    main:

    //
    // Map reads with STAR
    //
    def ch_star_out: Channel<StarAlignResult>
    if (use_sentieon_star) {
        ch_star_out = SENTIEON_STAR_ALIGN(ch_samples, index, gtf, star_ignore_sjdbgtf)
    } else if (use_parabricks_star) {
        ch_star_out = PARABRICKS_RNA_FQ2BAM(ch_samples, fasta, index, true, !skip_markduplicates)
    } else {
        ch_star_out = STAR_ALIGN(ch_samples, index, gtf, star_ignore_sjdbgtf)
    }

    // A run that produced no BAM drops out of the downstream channels.
    ch_star = ch_star_out.filter { r -> !r.raw_bams.isEmpty() }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_sorted: Channel<Bam> = BAM_SORT_STATS_SAMTOOLS(ch_star, fasta, fai)

    ch_star_sorted = ch_star.join(ch_sorted, by: 'id', remainder: true)
    ch_star_sorted.subscribe { r ->
        if( r.star == null || r.samtools == null ) {
            error "Sample '${r.id}' is missing its STAR sorted BAM result"
        }
    }
    def ch_results: Channel<StarAligned> = ch_star_sorted
        .filter { r -> r.star != null && r.samtools != null }
        .map { r -> r + record(aligner: 'star', percent_mapped: getStarPercentMapped(params, r.star.log_final)) }

    emit:
    ch_results
}
