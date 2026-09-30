nextflow.enable.types = true

//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'
include { ReadsInput; SalmonSubsampled } from '../../../modules/nf-core/types'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_samples: Channel<ReadsInput>
    ch_genome_fasta: Value<Path?> // decoys for the Salmon index, absent when the genome is not provided
    ch_transcript_fasta: Value<Path>
    ch_gtf: Value<Path>
    ch_index: Value<Path?> // prebuilt Salmon index, used when make_index is false
    make_index: Boolean // Whether to create salmon index before running salmon quant

    main:

    //
    // Create Salmon index if required
    //
    if (make_index) {
        ch_index_built = SALMON_INDEX(
            ch_transcript_fasta
                .map { transcript_fasta -> record(id: 'salmon_index', meta: [:], transcript_fasta: transcript_fasta) }
                .combine(genome_fasta: ch_genome_fasta)
        )
        ch_index_ref = ch_index_built.map { r -> r.index }
        ch_index_file = ch_index_ref
    }
    else {
        ch_index_ref = ch_index
        ch_index_file = ch_index.map { _index -> null }
    }

    //
    // Sub-sample FastQ files with fq
    //
    ch_subsampled = FQ_SUBSAMPLE(ch_samples)

    //
    // Pseudo-alignment with Salmon
    //
    ch_quant = SALMON_QUANT(ch_subsampled, ch_index_ref, ch_gtf, ch_transcript_fasta)

    ch_results = ch_subsampled.join(ch_quant, by: 'id')

    emit:
    samples     = ch_results
    index_built = ch_index_file // null unless the index was built here
}
