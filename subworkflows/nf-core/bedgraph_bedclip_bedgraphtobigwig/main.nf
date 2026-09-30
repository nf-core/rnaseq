//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { BigwigFiles           } from './types'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    //    path: chrom.sizes

    main:

    //
    // Clip bedGraph file
    //
    UCSC_BEDCLIP ( bedgraph, sizes )

    //
    // Convert bedGraph to bigWig
    //
    ch_bedgraph = UCSC_BEDCLIP.out.map { r -> [r.meta, r.bedgraph] }

    UCSC_BEDGRAPHTOBIGWIG ( ch_bedgraph, sizes )

    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG.out.map { r -> [r.meta, r.bigwig] }

    ch_results = ch_bigwig
        .join(ch_bedgraph)
        .map { meta, bigwig, clipped ->
            record(id: meta.id, bigwig: bigwig, bedgraph: clipped)
        }

    emit:
    bigwig   = ch_bigwig   // channel: [ val(meta), [ bigwig ] ]
    bedgraph = ch_bedgraph          // channel: [ val(meta), [ bedgraph ] ]
    results  = ch_results                       // channel: BigwigFiles
}
