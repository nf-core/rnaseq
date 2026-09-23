include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record Hisat2Logs {
    summary: Path
}

record Hisat2Aligned {
    id:       String
    meta:     Map
    aligner:  String
    orig_bam: Path
    unmapped: List<Path>?
    hisat2:   Hisat2Logs
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}

workflow FASTQ_ALIGN_HISAT2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/hisat2/index
    splicesites // channel: /path/to/genome.splicesites.txt
    ch_fasta_fai // channel: [meta, fasta, fai ]
    save_unaligned // val: boolean

    main:
    //
    // Map reads with HISAT2
    //
    HISAT2_ALIGN(reads, index, splicesites, save_unaligned)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(HISAT2_ALIGN.out.bam, ch_fasta_fai)

    ch_results = HISAT2_ALIGN.out.bam
        .join(HISAT2_ALIGN.out.summary)
        .join(HISAT2_ALIGN.out.fastq, remainder: true)
        .map { meta, orig_bam, summary, fastq ->
            record(
                id:       meta.id,
                meta:     meta,
                aligner:  'hisat2',
                orig_bam: orig_bam,
                unmapped: fastq ? [fastq].flatten() : null,
                hisat2:   record(summary: summary)
            )
        }
        .join(BAM_SORT_STATS_SAMTOOLS.out.results, by: 'id')
        .map { r -> r as Hisat2Aligned }

    emit:
    orig_bam = HISAT2_ALIGN.out.bam // channel: [ val(meta), bam   ]
    summary  = HISAT2_ALIGN.out.summary // channel: [ val(meta), log   ]
    fastq    = HISAT2_ALIGN.out.fastq // channel: [ val(meta), fastq ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam // channel: [ val(meta), [ bam ] ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index // channel: [ val(meta), [ index ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    results  = ch_results // channel: Hisat2Aligned
}
