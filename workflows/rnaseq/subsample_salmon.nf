//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'

include { getSalmonInferredStrandedness  } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_genome_fasta     // channel: /path/to/genome.fasta
    ch_transcript_fasta // channel: /path/to/transcript.fasta
    ch_gtf              // channel: /path/to/genome.gtf
    ch_index            // channel: /path/to/salmon/index/
    make_index          // boolean: Whether to create salmon index before running salmon quant

    main:
    //
    // Create Salmon index if required
    //
    def auto_reads = ch_reads
        .filter { meta, fastq -> meta.strandedness == 'auto' }
        .collect()

    def salmon_index = make_index && !auto_reads.isEmpty()
        ? SALMON_INDEX ( ch_genome_fasta, ch_transcript_fasta )
        : ch_index

    ch_trim_reads = ch_reads.map { meta, fastq ->

        if( meta.strandedness != 'auto' )
            return [ meta, fastq ]

        //
        // Sub-sample FastQ files with fq
        //
        def subsamples = FQ_SUBSAMPLE ( meta, fastq )

        //
        // Pseudo-alignment with Salmon
        //
        def lib_type = 'A'
        def alignment_mode = false
        def (results, json_info) = SALMON_QUANT ( subsamples, salmon_index, ch_gtf, ch_transcript_fasta, alignment_mode, lib_type )

        [ meta + [ strandedness: getSalmonInferredStrandedness(json_info) ], subsamples ]
    }

    emit:
    ch_trim_reads

}
