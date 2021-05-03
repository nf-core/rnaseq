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

params.fasta        = WorkflowRnaseq.getGenomeAttribute(params, 'fasta')
params.gtf          = WorkflowRnaseq.getGenomeAttribute(params, 'gtf')
params.gff          = WorkflowRnaseq.getGenomeAttribute(params, 'gff')
params.gene_bed     = WorkflowRnaseq.getGenomeAttribute(params, 'bed12')
params.star_index   = WorkflowRnaseq.getGenomeAttribute(params, 'star')
params.hisat2_index = WorkflowRnaseq.getGenomeAttribute(params, 'hisat2')
params.rsem_index   = WorkflowRnaseq.getGenomeAttribute(params, 'rsem')
params.salmon_index = WorkflowRnaseq.getGenomeAttribute(params, 'salmon')

////////////////////////////////////////////////////
/* --    INITIALISE & CHECK MAIN PARAMETERS    -- */
////////////////////////////////////////////////////
    
WorkflowMain.initialise(workflow, params, log)

////////////////////////////////////////////////////
/* --        NAMED PIPELINE WORKFLOW           -- */
////////////////////////////////////////////////////

workflow NFCORE_RNASEQ {
    
    /*
     * WORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
     */
    if (params.public_data_ids) {
        include { SRA_DOWNLOAD } from './workflows/sra_download'
        SRA_DOWNLOAD ()
    
    /*
     * WORKFLOW: Run main nf-core/rnaseq analysis pipeline
     */
    } else {
        include { RNASEQ } from './workflows/rnaseq'
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
