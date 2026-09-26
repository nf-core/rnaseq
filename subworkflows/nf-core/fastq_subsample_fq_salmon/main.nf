//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_genome_fasta     // channel: /path/to/genome.fasta
    ch_transcript_fasta // channel: /path/to/transcript.fasta
    ch_gtf              // channel: /path/to/genome.gtf
    ch_index            // channel: /path/to/salmon/index/
    make_index          // boolean: Whether to create salmon index before running salmon quant

    main:

    ch_gtf_transcript_fasta = ch_gtf.combine(ch_transcript_fasta)

    //
    // Create Salmon index if required
    //
    ch_index_built = channel.empty()
    if (make_index) {
        // genome_fasta may be an empty list (no decoys); wrap before combine() so an
        // empty list contributes a position instead of being flattened away
        ch_transcript_fasta
            .combine(ch_genome_fasta.map { genome_fasta -> [genome_fasta] })
            .map { items -> [ [:], items[0], items[1] ] }
            .set { ch_index_input }

        ch_index = SALMON_INDEX ( ch_index_input ).index
        ch_index_built = ch_index
    }
    else {
        ch_index = ch_index.map { index -> [ [:], index ] }
    }

    //
    // Sub-sample FastQ files with fq
    //
    FQ_SUBSAMPLE ( ch_reads )

    //
    // Pseudo-alignment with Salmon
    //
    SALMON_QUANT ( FQ_SUBSAMPLE.out.fastq, ch_index.combine(ch_gtf_transcript_fasta).first() )

    emit:
    index             = ch_index                           // channel: [ val(meta), index ]
    index_built       = ch_index_built.map { _meta, index -> index } // channel: path(salmon/index/), only set when built here

    reads             = FQ_SUBSAMPLE.out.fastq             // channel: [ val(meta), fastq ]

    results           = SALMON_QUANT.out.results           // channel: [ val(meta), results_dir ]
    json_info         = SALMON_QUANT.out.json_info         // channel: [ val(meta), json_info
    lib_format_counts = SALMON_QUANT.out.lib_format_counts // channel: [ val(meta), json_info
}
