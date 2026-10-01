nextflow.enable.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'
include { BamInput; SamtoolsStats; MarkdupBam } from '../../../modules/nf-core/types'
include { PicardMarkduplicatesResult } from '../../../modules/nf-core/picard/markduplicates/main'
include { SamtoolsIndexResult } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam: Channel<BamInput>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    run_stats: Boolean
    index_args: String // samtools index options

    main:
    def ch_markdup: Channel<PicardMarkduplicatesResult> = PICARD_MARKDUPLICATES(ch_bam, ch_fasta, ch_fai)

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    ch_marked = ch_markdup.map { r -> record(id: r.id, meta: r.meta, bam: r.bam ?: r.cram) }

    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_marked.map { r -> r + record(args: index_args) })

    def ch_stats: Channel<SamtoolsStats> = channel.empty()
    if (run_stats) {
        ch_stats = BAM_STATS_SAMTOOLS(ch_marked.join(ch_index, by: 'id'), ch_fasta, ch_fai)
    }

    ch_markdup_indexed = ch_markdup.join(ch_index, by: 'id', remainder: true)
    ch_markdup_indexed.subscribe { r ->
        if( r.metrics == null || r.bai == null ) {
            error "Sample '${r.id}' is missing its samtools index result"
        }
    }
    ch_results = ch_markdup_indexed
        .filter { r -> r.metrics != null && r.bai != null }
        // remainder keeps samples without samtools stats when the caller disabled them.
        .join(ch_stats, by: 'id', remainder: true)

    emit:
    ch_results
}
