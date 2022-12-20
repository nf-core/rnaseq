//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_index            // channel: /path/to/salmon/index/
    ch_transcript_fasta // channel: /path/to/transcript.fasta
    ch_gtf              // channel: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Sub-sample FastQ files with fq
    //
    FQ_SUBSAMPLE ( ch_reads )
    ch_versions = ch_versions.mix(FQ_SUBSAMPLE.out.versions.first())

    //
    // Pseudo-alignment with Salmon
    //
    def lib_type = 'A'
    def alignment_mode = false
    SALMON_QUANT ( FQ_SUBSAMPLE.out.fastq, ch_index, ch_gtf, ch_transcript_fasta, alignment_mode, lib_type )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
    reads     = FQ_SUBSAMPLE.out.fastq     // channel: [ val(meta), fastq ]

    results   = SALMON_QUANT.out.results   // channel: [ val(meta), results_dir ]
    json_info = SALMON_QUANT.out.json_info // channel: [ val(meta), json_info

    versions  = ch_versions                // channel: [ versions.yml ]
}
