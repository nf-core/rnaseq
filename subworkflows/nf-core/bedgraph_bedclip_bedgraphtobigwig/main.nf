//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

record BigwigFiles {
    id:       String
    bigwig:   Path
    bedgraph: Path
}

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
    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.bedgraph, sizes )

    ch_results = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
        .join(UCSC_BEDCLIP.out.bedgraph)
        .map { meta, bigwig, clipped ->
            record(id: meta.id, bigwig: bigwig, bedgraph: clipped) as BigwigFiles
        }

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    results  = ch_results                       // channel: BigwigFiles
}
