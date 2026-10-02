nextflow.enable.types = true

//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { BedgraphInput; BigwigFiles } from '../../../modules/nf-core/types'
include { UcscBedclipResult } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UcscBedgraphtobigwigResult } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    ch_bedgraph: Channel<BedgraphInput>
    sizes: Value<Path> // chrom.sizes

    main:

    //
    // Clip bedGraph file
    //
    def ch_clipped: Channel<UcscBedclipResult> = UCSC_BEDCLIP(ch_bedgraph, sizes)

    //
    // Convert bedGraph to bigWig
    //
    def ch_bigwig: Channel<UcscBedgraphtobigwigResult> = UCSC_BEDGRAPHTOBIGWIG(ch_clipped, sizes)

    ch_results = ch_clipped.join(ch_bigwig, by: 'id')

    emit:
    ch_results
}
