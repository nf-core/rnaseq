include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'

include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index'

include { SAMTOOLS_STATS          } from '../../../modules/nf-core/samtools/stats'
include { SAMTOOLS_IDXSTATS       } from '../../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT       } from '../../../modules/nf-core/samtools/flagstat'

include { UMITOOLS_DEDUP          } from '../../../modules/nf-core/umitools/dedup'

workflow FASTQ_ALIGN_HISAT2 {

    take:
    reads       // channel: [ meta: Map, fastq: List<Path> ]
    index       // path: /path/to/hisat2/index
    splicesites // path: /path/to/genome.splicesites.txt
    fasta       // path

    main:

    bam = reads.map { meta, fastq ->
        //
        // Map reads with HISAT2
        //
        def hisat2 = HISAT2_ALIGN ( meta, fastq, index, splicesites )

        //
        // Sort, index BAM file and run samtools stats, flagstat and idxstats
        //
        def bam_sort = SAMTOOLS_SORT ( meta, hisat2.bam, fasta )
        def bam_index = SAMTOOLS_INDEX ( meta, bam_sort.bam )
        def bam_bai_csi = bam_index.bai ?: bam_index.csi
        def bam_stats = [
            SAMTOOLS_STATS ( meta, bam_index.bam, bam_bai_csi, fasta ),
            SAMTOOLS_FLAGSTAT ( meta, bam_index.bam, bam_bai_csi ),
            SAMTOOLS_IDXSTATS ( meta, bam_index.bam, bam_bai_csi ),
        ]

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            // Deduplicate genome BAM file before downstream analysis
            def umi_dedup = UMITOOLS_DEDUP ( meta, bam_index.bam, bam_index.bai, params.umitools_dedup_stats )
            def umi_index = SAMTOOLS_INDEX ( meta, umi_dedup.bam )
            def umi_bai_csi = umi_index.bai ?: umi_index.csi
            def umi_stats = [
                SAMTOOLS_STATS ( meta, umi_index.bam, umi_bai_csi, fasta ),
                SAMTOOLS_FLAGSTAT ( meta, umi_index.bam, umi_bai_csi ),
                SAMTOOLS_IDXSTATS ( meta, umi_index.bam, umi_bai_csi ),
            ]

            // TODO: add umi results to record type
            [ meta, umi_dedup, umi_index, umi_stats ]
        }

        // TODO: Hisat2Sample record type
        [ meta, hisat2, bam_stats ]
    }

    emit:
    bam

}

