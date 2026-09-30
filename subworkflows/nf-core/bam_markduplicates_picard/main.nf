nextflow.enable.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'
include { BamToMarkdup; MarkdupBam } from './types'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam: Channel<BamToMarkdup>
    ch_fasta_fai: Value<Tuple<Map, Path?, Path?>>
    run_stats: Boolean

    main:
    ch_markdup = PICARD_MARKDUPLICATES(ch_bam, ch_fasta_fai)

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    ch_marked = ch_markdup.map { r -> record(id: r.id, meta: r.meta, bam: r.bam ?: r.cram) }

    ch_index = SAMTOOLS_INDEX(ch_marked)

    ch_stats = channel.empty()
    if (run_stats) {
        ch_stats = BAM_STATS_SAMTOOLS(ch_marked.join(ch_index, by: 'id'), ch_fasta_fai)
    }

    // remainder keeps samples without samtools stats when the caller disabled them.
    ch_results = ch_markdup
        .join(ch_index, by: 'id')
        .join(ch_stats, by: 'id', remainder: true)

    emit:
    ch_results
}
