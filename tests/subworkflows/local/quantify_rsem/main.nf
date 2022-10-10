#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()

include { QUANTIFY_RSEM  } from '../../../../subworkflows/local/quantify_rsem'  addParams(
    calculateexpression_options: rsem_calculateexpression_options,
    samtools_sort_options: samtools_sort_genome_options,
    samtools_index_options: samtools_index_genome_options,
    samtools_stats_options: samtools_index_genome_options,
    merge_counts_options: modules['rsem_merge_counts']
)

workflow test_quantify_rsem {
    input = [ [ id:'test', single_end:true ],
             file(params.test_data['sarscov2']['illumina']['bam']['test_single_end_bam'], checkIfExists: true)]

    rsem_index = file('https://github.com/nf-core/test-datasets/raw/rnaseq/reference/rsem.tar.gz', checkIfExists: true)

    QUANTIFY_RSEM (
        input,
        rsem_index
    )
}
