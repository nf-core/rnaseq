//
// Clip over-running ends from bedGraph file and convert to bigWig
//

params.bedclip_options          = [:]
params.bedgraphtobigwig_options = [:]

include { UCSC_BEDCLIP          } from '../../modules/nf-core/modules/ucsc/bedclip/main'          addParams( options: params.bedclip_options          )
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main' addParams( options: params.bedgraphtobigwig_options )

workflow BEDGRAPH_TO_BIGWIG {
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

    emit:
    bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph     = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    ucsc_version = UCSC_BEDCLIP.out.version         //    path: versions.yml
}
