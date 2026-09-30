//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'
include { MarkdupBam            } from './types'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_reads // channel: [ val(meta), path(reads) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai)]
    run_stats // boolean: run samtools stats, flagstat and idxstats

    main:
    PICARD_MARKDUPLICATES(ch_reads, ch_fasta_fai)

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    ch_markdup = PICARD_MARKDUPLICATES.out.map { r -> [r.meta, r.bam ?: r.cram] }

    SAMTOOLS_INDEX(ch_markdup)

    ch_reads_index = ch_markdup.join(SAMTOOLS_INDEX.out, by: [0])

    ch_stats = channel.empty()
    ch_flagstat = channel.empty()
    ch_idxstats = channel.empty()
    ch_samtools = channel.empty()
    if (run_stats) {
        BAM_STATS_SAMTOOLS(ch_reads_index, ch_fasta_fai)
        ch_stats = BAM_STATS_SAMTOOLS.out.stats
        ch_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
        ch_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
        ch_samtools = BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }
    }

    // samtools is a remainder join because the stats can be disabled by the caller.
    ch_results = PICARD_MARKDUPLICATES.out.map { r -> [r.id, r] }
        .join(SAMTOOLS_INDEX.out.map { meta, index -> [meta.id, index] }, by: [0])
        .join(ch_samtools, by: [0], remainder: true)
        .map { id, markdup, index, samtools ->
            record(
                id: id,
                bam: markdup.bam,
                cram: markdup.cram,
                bai: index,
                metrics: markdup.metrics,
                samtools: samtools ? record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats) : null
            )
        }

    ch_metrics = PICARD_MARKDUPLICATES.out.map { r -> [r.meta, r.metrics] }

    ch_per_sample_mqc_bundle = ch_stats
        .join(ch_flagstat,   remainder: true)
        .join(ch_idxstats,   remainder: true)
        .join(ch_metrics,    remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    emit:
    bam                   = PICARD_MARKDUPLICATES.out.filter { r -> r.bam }.map { r -> [r.meta, r.bam] } // channel: [ val(meta), path(bam) ]
    cram                  = PICARD_MARKDUPLICATES.out.filter { r -> r.cram }.map { r -> [r.meta, r.cram] } // channel: [ val(meta), path(cram) ]
    metrics               = ch_metrics // channel: [ val(meta), path(metrics) ]
    index                 = SAMTOOLS_INDEX.out // channel: [ val(meta), path(index) ]
    stats                 = ch_stats // channel: [ val(meta), path(stats) ]
    flagstat              = ch_flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats              = ch_idxstats // channel: [ val(meta), path(idxstats) ]
    per_sample_mqc_bundle = ch_per_sample_mqc_bundle // channel: [ val(meta), list(files) ]
    results               = ch_results // channel: MarkdupBam
}
