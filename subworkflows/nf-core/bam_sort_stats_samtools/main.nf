nextflow.preview.types = true

//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

record SamtoolsResult {
    meta:     Map
    bam:      Path
    bai:      Path?
    csi:      Path?
    stats:    Path
    flagstat: Path
    idxstats: Path
}

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = channel.empty()

    SAMTOOLS_SORT ( ch_bam, ch_fasta, '' )

    // Extract bam from SamtoolsSortResult record
    ch_sorted_bam = SAMTOOLS_SORT.out.map { r -> [r.meta, r.bam] }

    SAMTOOLS_INDEX ( ch_sorted_bam )

    ch_index_bai = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.bai] }
    ch_index_csi = SAMTOOLS_INDEX.out.map { r -> [r.meta, r.csi] }

    ch_sorted_bam
        .join(ch_index_bai, by: [0], remainder: true)
        .join(ch_index_csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    // Aggregate outputs from multiple processes into a single record channel
    // using join (by meta key) + map to construct the record.
    emit:
    result = ch_sorted_bam
        .join(ch_index_bai, by: [0], remainder: true)
        .join(ch_index_csi, by: [0], remainder: true)
        .join(BAM_STATS_SAMTOOLS.out.stats,    by: [0])
        .join(BAM_STATS_SAMTOOLS.out.flagstat,  by: [0])
        .join(BAM_STATS_SAMTOOLS.out.idxstats,  by: [0])
        .map { meta, bam, bai, csi, stats, flagstat, idxstats ->
            record(
                meta: meta,
                bam: bam, bai: bai, csi: csi,
                stats: stats, flagstat: flagstat, idxstats: idxstats
            )
        }
    versions = ch_versions                     // channel: [ versions.yml ]
}
