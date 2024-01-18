//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    //    path: chrom.sizes

    main:

    ( bedgraph & sizes )
        | UCSC_BEDCLIP
        | { out -> UCSC_BEDGRAPHTOBIGWIG ( out.bedgraph, sizes ) }

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}
