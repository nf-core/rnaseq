//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {

    take:
    reads   // channel: [ meta: Map, input: Path ]
    fasta   // file
    fai     // file

    main:

    markdup = reads.map { meta, input ->
        def out = PICARD_MARKDUPLICATES ( meta, input, fasta, fai )
        [ meta, out ]
    }

    markdup_indexed = markdup
        .flatMap { meta, out ->
            def result = []
            if( out.bam )
                result << [ meta, out.bam ]
            if( out.cram )
                result << [ meta, out.cram ]
            result
        }
        .map { meta, bam_cram ->
            def index = SAMTOOLS_INDEX ( meta, bam_cram )
            [ meta, bam_cram, index.bai ?: index.crai ?: index.csi ]
        }

    indexed = markdup_indexed.map { meta, _bam_cram, index -> [ meta, index ] }

    stats = BAM_STATS_SAMTOOLS (
        markdup_indexed,
        fasta
    )

    emit:
    markdup     // channel: [ meta: Map, bam: Path, cram: Path, metrics: Path ]
    indexed     // channel: [ meta: Map, bai: Path, crai: Path, csi: Path ]
    stats       // channel: [ meta: Map, stats: Path, flagstat: Path, idxstats: Path ]

}
