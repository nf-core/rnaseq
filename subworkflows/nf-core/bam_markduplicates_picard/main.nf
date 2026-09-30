nextflow.enable.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'
include { PicardMarkduplicatesResult; SamtoolsStats; SamtoolsIndexResult; Bam; MarkdupBam } from '../../../modules/nf-core/types'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam: Channel<Bam>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    run_stats: Boolean

    main:
    def ch_markdup: Channel<PicardMarkduplicatesResult> = PICARD_MARKDUPLICATES(ch_bam, ch_fasta, ch_fai)

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    ch_marked = ch_markdup.map { r -> record(id: r.id, meta: r.meta, bam: r.bam ?: r.cram) }

    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_marked)

    def ch_stats: Channel<SamtoolsStats> = channel.empty()
    if (run_stats) {
        ch_stats = BAM_STATS_SAMTOOLS(ch_marked.join(ch_index, by: 'id'), ch_fasta, ch_fai)
    }

    // remainder keeps samples without samtools stats when the caller disabled them.
    def ch_results: Channel<MarkdupBam> = ch_markdup
        .join(ch_index, by: 'id')
        .join(ch_stats, by: 'id', remainder: true)

    emit:
    ch_results
}
