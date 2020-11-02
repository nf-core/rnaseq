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

def json_schema = "$baseDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

// Force print these hidden parameters in the JSON Schema
def force_params = [
    'max_memory', 'max_cpus', 'max_time',
    'config_profile_description', 'config_profile_contact', 'config_profile_url'
]
summary = Schema.params_summary_map(workflow, params, json_schema, force_params)
log.info  Schema.params_summary_log(workflow, params, json_schema, force_params)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// workflow {
//     if (params.public_data_ids) {
//         /*
//          * SUBWORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
//          */
//         include { SRA_DOWNLOAD } from './sra_download'
//         SRA_DOWNLOAD ()
//     } else {
//         /*
//          * SUBWORKFLOW: Run main nf-core/rnaseq analysis pipeline
//          */
//         include { RNASEQ } from './rnaseq' addParams( summary: summary )
//         RNASEQ ()
//     }
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
