//
// Subworkflow that uses the nf-schema plugin to validate parameters and render the parameter summary
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    input_workflow      // workflow: the workflow object used by nf-schema to get metadata from the workflow
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   path to the parameters JSON schema.
                        //           this has to be the same as the schema given to `validation.parametersSchema`
                        //           when this input is empty it will automatically use the configured schema or
                        //           "${projectDir}/nextflow_schema.json" as default. This input should not be empty
                        //           for meta pipelines

    main:

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //
    if(parameters_schema) {
        log.info paramsSummaryLog(input_workflow, parameters_schema:parameters_schema)
    } else {
        log.info paramsSummaryLog(input_workflow)
    }

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(validate_params) {
        if(parameters_schema) {
            validateParameters(parameters_schema:parameters_schema)
        } else {
            validateParameters()
        }
    }

    emit:
    dummy_emit = true
}

