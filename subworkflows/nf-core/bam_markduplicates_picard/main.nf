nextflow.preview.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

record MarkDupResult {
    meta:     Map
    bam:      Path
    bai:      Path
    metrics:  Path
    stats:    Path
    flagstat: Path
    idxstats: Path
}

workflow BAM_MARKDUPLICATES_PICARD {

    take:
    ch_reads   // channel: [ val(meta), path(reads) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai   // channel: [ val(meta), path(fai) ]

    main:

    ch_versions = channel.empty()

    PICARD_MARKDUPLICATES ( ch_reads, ch_fasta, ch_fai )

    // Extract bam/cram from the record for downstream indexing
    ch_picard_bam  = PICARD_MARKDUPLICATES.out.map { r -> [r.meta, r.bam] }.filter { it[1] }
    ch_picard_cram = PICARD_MARKDUPLICATES.out.map { r -> [r.meta, r.cram] }.filter { it[1] }
    ch_markdup = ch_picard_bam.mix(ch_picard_cram)

    SAMTOOLS_INDEX ( ch_markdup )

    // Extract index fields from SamtoolsIndexResult record
    ch_index_bai  = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.bai] }
    ch_index_crai = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.crai] }
    ch_index_csi  = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.csi] }

    ch_reads_index = ch_markdup
        .join(ch_index_bai,  by: [0], remainder: true)
        .join(ch_index_crai, by: [0], remainder: true)
        .join(ch_index_csi,  by: [0], remainder: true)
        .map{meta, reads, bai, crai, csi ->
            if (bai) [ meta, reads, bai ]
            else if (crai) [ meta, reads, crai ]
            else [ meta, reads, csi ]
        }

    BAM_STATS_SAMTOOLS ( ch_reads_index, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    ch_picard_metrics = PICARD_MARKDUPLICATES.out.map { r -> [r.meta, r.metrics] }

    emit:
    result = ch_picard_bam
        .join(ch_index_bai, by: [0], remainder: true)
        .join(ch_picard_metrics,              by: [0])
        .join(BAM_STATS_SAMTOOLS.out.stats,    by: [0])
        .join(BAM_STATS_SAMTOOLS.out.flagstat,  by: [0])
        .join(BAM_STATS_SAMTOOLS.out.idxstats,  by: [0])
        .map { meta, bam, bai, metrics, stats, flagstat, idxstats ->
            record(
                meta: meta,
                bam: bam, bai: bai, metrics: metrics,
                stats: stats, flagstat: flagstat, idxstats: idxstats
            )
        }
    bam      = ch_picard_bam
    csi      = ch_index_csi
    metrics  = ch_picard_metrics
    stats    = BAM_STATS_SAMTOOLS.out.stats
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats
    versions = ch_versions
}
