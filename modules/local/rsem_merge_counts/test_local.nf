#! /usr/bin/env nextflow

// /home/eladh/code/development/NF-core/nf-core-rnaseq_3.21.0/3_21_0/modules/local/rsem_merge_counts/genes
// /home/eladh/code/development/NF-core/nf-core-rnaseq_3.21.0/3_21_0/modules/local/rsem_merge_counts/isoforms

nextflow.enable.dsl=2
include {RSEM_MERGE_COUNTS} from './main.nf'



workflow  {

    RSEM_MERGE_COUNTS(Channel.fromPath(params.genes), Channel.fromPath(params.isoforms))


}
