//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record MarkdupBam {
    id:       String
    bam:      Path?
    cram:     Path?
    bai:      Path
    metrics:  Path
    samtools: SamtoolsStatsFiles?
}

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_reads // channel: [ val(meta), path(reads) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai)]

    main:
    PICARD_MARKDUPLICATES(ch_reads, ch_fasta_fai)

    ch_markdup = PICARD_MARKDUPLICATES.out.bam.mix(PICARD_MARKDUPLICATES.out.cram)

    SAMTOOLS_INDEX(ch_markdup)

    ch_reads_index = ch_markdup.join(SAMTOOLS_INDEX.out.index, by: [0])

    BAM_STATS_SAMTOOLS(ch_reads_index, ch_fasta_fai)

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    // BAM_STATS_SAMTOOLS may be disabled by the caller via ext.when, so samtools is a remainder join.
    ch_results = PICARD_MARKDUPLICATES.out.bam.map { meta, bam -> [meta.id, bam, null] }
        .mix(PICARD_MARKDUPLICATES.out.cram.map { meta, cram -> [meta.id, null, cram] })
        .join(SAMTOOLS_INDEX.out.index.map { meta, index -> [meta.id, index] }, by: [0])
        .join(PICARD_MARKDUPLICATES.out.metrics.map { meta, metrics -> [meta.id, metrics] }, by: [0])
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0], remainder: true)
        .map { id, bam, cram, index, metrics, samtools ->
            record(
                id: id,
                bam: bam,
                cram: cram,
                bai: index,
                metrics: metrics,
                samtools: samtools ? record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats) : null
            ) as MarkdupBam
        }

    ch_per_sample_mqc_bundle = BAM_STATS_SAMTOOLS.out.stats
        .join(BAM_STATS_SAMTOOLS.out.flagstat,   remainder: true)
        .join(BAM_STATS_SAMTOOLS.out.idxstats,   remainder: true)
        .join(PICARD_MARKDUPLICATES.out.metrics, remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    emit:
    bam                   = PICARD_MARKDUPLICATES.out.bam // channel: [ val(meta), path(bam) ]
    cram                  = PICARD_MARKDUPLICATES.out.cram // channel: [ val(meta), path(cram) ]
    metrics               = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(metrics) ]
    index                 = SAMTOOLS_INDEX.out.index // channel: [ val(meta), path(index) ]
    stats                 = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat              = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats              = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    per_sample_mqc_bundle = ch_per_sample_mqc_bundle // channel: [ val(meta), list(files) ]
    results               = ch_results // channel: MarkdupBam
}
