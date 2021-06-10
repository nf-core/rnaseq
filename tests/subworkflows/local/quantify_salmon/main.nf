#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()

include { SALMON_INDEX } from '../../../../modules/nf-core/software/salmon/index/main.nf' addParams( options: [:] )
include { QUANTIFY_SALMON } from '../../../../subworkflows/local/quantify_salmon' addParams(
    genome_options: '',
    tximport_options: modules['salmon_tximport'],
    salmon_quant_options: salmon_quant_options,
    merge_counts_options: modules['salmon_merge_counts']
)

workflow test_quantify_salmon {
    input = [ [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]
    genome_fasta     = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    transcript_fasta = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    gtf              = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    QUANTIFY_SALMON ( input, SALMON_INDEX.out.index, transcript_fasta, gtf, false )
}
