nextflow.enable.types = true

//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { BedgraphInput; BigwigFiles } from '../../../modules/nf-core/types'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    ch_bedgraph: Channel<BedgraphInput>
    sizes: Value<Path> // chrom.sizes

    main:

    //
    // Clip bedGraph file
    //
    ch_clipped = UCSC_BEDCLIP(ch_bedgraph, sizes)

    //
    // Convert bedGraph to bigWig
    //
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG(ch_clipped, sizes)

    ch_results = ch_clipped.join(ch_bigwig, by: 'id')

    emit:
    ch_results
}
