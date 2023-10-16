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

    ch_versions = Channel.empty()

    //
    // Create Salmon index if required
    //
    if (make_index) {
        SALMON_INDEX.config.ext.args   = params.gencode ? '--gencode' : ''
        SALMON_INDEX.config.publishDir = [
            path: "${params.outdir}/genome/index",
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
            enabled: params.save_reference
        ]
        ch_index = SALMON_INDEX ( ch_genome_fasta, ch_transcript_fasta ).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    }

    //
    // Sub-sample FastQ files with fq
    //
    FQ_SUBSAMPLE.config.ext.args   = params.subsample_fq_args
    FQ_SUBSAMPLE.config.ext.prefix = { "${meta.id}.subsampled" }
    FQ_SUBSAMPLE.config.publishDir = [
        path: "${params.outdir}/sample_fastq/fastq",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        enabled: false
    ]
    FQ_SUBSAMPLE ( ch_reads )
    ch_versions = ch_versions.mix(FQ_SUBSAMPLE.out.versions.first())

    //
    // Pseudo-alignment with Salmon
    //
    def lib_type = 'A'
    def alignment_mode = false
    SALMON_QUANT.config.ext.args   = params.subsample_salmon_args
    SALMON_QUANT.config.publishDir = [
        path: "${params.outdir}/sample_fastq/salmon",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') || filename.endsWith('_meta_info.json') ? null : filename },
        enabled: false
    ]
    SALMON_QUANT ( FQ_SUBSAMPLE.out.fastq, ch_index, ch_gtf, ch_transcript_fasta, alignment_mode, lib_type )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
    index     = ch_index                   // channel: [ index ]

    reads     = FQ_SUBSAMPLE.out.fastq     // channel: [ val(meta), fastq ]

    results   = SALMON_QUANT.out.results   // channel: [ val(meta), results_dir ]
    json_info = SALMON_QUANT.out.json_info // channel: [ val(meta), json_info

    versions  = ch_versions                // channel: [ versions.yml ]
}
