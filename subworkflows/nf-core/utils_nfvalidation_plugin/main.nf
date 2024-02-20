//
// Subworkflow that uses the nf-validation plugin to render help text and parameter summary
//

/*
========================================================================================
    IMPORT NF-VALIDATION PLUGIN
========================================================================================
*/

include { paramsHelp         } from 'plugin/nf-validation'
include { paramsSummaryLog   } from 'plugin/nf-validation'
include { validateParameters } from 'plugin/nf-validation'

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFVALIDATION_PLUGIN {

    take:
    print_help       // boolean: print help
    workflow_command //  string: default commmand used to run pipeline
    pre_help_text    //  string: string to be printed before help text and summary log
    post_help_text   //  string: string to be printed after help text and summary log
    validate_params  // boolean: validate parameters
    schema_filename  //    path: JSON schema file, null to use default value

    main:

    log.debug "Using schema file: ${schema_filename}"

    // Default values for strings
    pre_help_text    = pre_help_text    ?: ''
    post_help_text   = post_help_text   ?: ''
    workflow_command = workflow_command ?: ''

    //
    // Print help message if needed
    //
    if (print_help) {
        log.info pre_help_text + paramsHelp(workflow_command, parameters_schema: schema_filename) + post_help_text
        System.exit(0)
    }

    //
    // Print parameter summary to stdout
    //
    log.info pre_help_text + paramsSummaryLog(workflow, parameters_schema: schema_filename) + post_help_text

    //
    // Validate parameters relative to the parameter JSON schema
    //
    if (validate_params){
        validateParameters(parameters_schema: schema_filename)
    }

    emit:
    dummy_emit = true
}
