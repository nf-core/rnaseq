#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnaseq
========================================================================================
 nf-core/rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.fasta        = Workflow.getGenomeAttribute(params, 'fasta')
params.gtf          = Workflow.getGenomeAttribute(params, 'gtf')
params.gff          = Workflow.getGenomeAttribute(params, 'gff')
params.gene_bed     = Workflow.getGenomeAttribute(params, 'bed12')
params.star_index   = Workflow.getGenomeAttribute(params, 'star')
params.hisat2_index = Workflow.getGenomeAttribute(params, 'hisat2')
params.rsem_index   = Workflow.getGenomeAttribute(params, 'rsem')
params.salmon_index = Workflow.getGenomeAttribute(params, 'salmon')

////////////////////////////////////////////////////
/* --          VALIDATE PARAMETERS             -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --        NAMED PIPELINE WORKFLOW           -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params, json_schema)

workflow NFCORE_RNASEQ {
    
    /*
     * WORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
     */
    if (params.public_data_ids) {
        include { SRA_DOWNLOAD } from './workflows/sra_download' addParams( summary_params: summary_params )
        SRA_DOWNLOAD ()
    
    /*
     * WORKFLOW: Run main nf-core/rnaseq analysis pipeline
     */
    } else {
        include { RNASEQ } from './workflows/rnaseq' addParams( summary_params: summary_params )
        RNASEQ ()
    }
}

////////////////////////////////////////////////////
/* --           RUN ALL WORKFLOWS              -- */
////////////////////////////////////////////////////

/*
 * WORKFLOW: Execute a single named workflow for the pipeline
 * See: https://github.com/nf-core/rnaseq/issues/619
 */
workflow {
    NFCORE_RNASEQ ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
