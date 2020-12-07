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
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.fasta        = Checks.get_genome_attribute(params, 'fasta')
params.gtf          = Checks.get_genome_attribute(params, 'gtf')
params.gff          = Checks.get_genome_attribute(params, 'gff')
params.gene_bed     = Checks.get_genome_attribute(params, 'bed12')
params.star_index   = Checks.get_genome_attribute(params, 'star')
params.hisat2_index = Checks.get_genome_attribute(params, 'hisat2')
params.rsem_index   = Checks.get_genome_attribute(params, 'rsem')
params.salmon_index = Checks.get_genome_attribute(params, 'salmon')

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    if (params.public_data_ids) {
        /*
         * SUBWORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
         */
        include { SRA_DOWNLOAD } from './sra_download' addParams( summary_params: summary_params )
        SRA_DOWNLOAD ()
    } else {
        /*
         * SUBWORKFLOW: Run main nf-core/rnaseq analysis pipeline
         */
        include { RNASEQ } from './rnaseq' addParams( summary_params: summary_params )
        RNASEQ ()
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////