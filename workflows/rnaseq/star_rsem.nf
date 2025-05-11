//
// Alignment with STAR and gene/transcript quantification with Salmon
//

include { RSEM_CALCULATEEXPRESSION } from '../../modules/nf-core/rsem/calculateexpression'
include { RSEM_MERGE_COUNTS        } from '../../modules/local/rsem_merge_counts'

include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index'

include { SAMTOOLS_STATS          } from '../../modules/nf-core/samtools/stats'
include { SAMTOOLS_IDXSTATS       } from '../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT       } from '../../modules/nf-core/samtools/flagstat'


workflow STAR_RSEM {
    take:
    reads   // channel: [ val(meta), [ reads ] ]
    index   // path: /path/to/rsem/index/
    fasta   // path: /path/to/fasta

    main:

    bam = reads.map { meta, fastq ->
        //
        // Quantify reads with RSEM
        //
        def rsem = RSEM_CALCULATEEXPRESSION ( meta, fastq, index )

        //
        // Sort, index BAM file and run samtools stats, flagstat and idxstats
        //
        def bam_sort = SAMTOOLS_SORT ( meta, rsem.bam_star, fasta )
        def bam_index = SAMTOOLS_INDEX ( meta, bam_sort.bam )
        def bam_bai_csi = bam_index.bai ?: bam_index.csi
        def bam_stats = [
            SAMTOOLS_STATS ( meta, bam_index.bam, bam_bai_csi, fasta ),
            SAMTOOLS_FLAGSTAT ( meta, bam_index.bam, bam_bai_csi ),
            SAMTOOLS_IDXSTATS ( meta, bam_index.bam, bam_bai_csi ),
        ]

        // TODO: StarRsemSample record type
        [ meta, rsem, bam_stats ]
    }

    // TODO: StarRsemCounts record type
    //
    // Merge counts across samples
    //
    counts = RSEM_MERGE_COUNTS (
        bam.collect { sample -> sample.counts_gene },
        bam.collect { sample -> sample.counts_transcript }
    )

    emit:
    bam
    counts

}
