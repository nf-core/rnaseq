//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record SortedBam {
    id:       String
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]

    main:
    SAMTOOLS_SORT(ch_bam, ch_fasta_fai, '')

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.index, by: [0])
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    ch_results = ch_bam_bai
        .map { meta, bam, bai -> [meta.id, bam, bai] }
        .join(BAM_STATS_SAMTOOLS.out.results.map { r -> [r.id, r] }, by: [0])
        .map { id, bam, bai, samtools ->
            record(
                id: id,
                bam: bam,
                bai: bai,
                samtools: record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats)
            ) as SortedBam
        }

    emit:
    bam      = SAMTOOLS_SORT.out.bam // channel: [ val(meta), [ bam ] ]
    index    = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    results  = ch_results // channel: SortedBam
}
