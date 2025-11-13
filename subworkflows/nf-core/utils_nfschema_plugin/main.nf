//
// Subworkflow that uses the nf-schema plugin to validate parameters and render the parameter summary
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsHelp         } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    input_workflow      // workflow: the workflow object used by nf-schema to get metadata from the workflow
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   path to the parameters JSON schema.
                        //           this has to be the same as the schema given to `validation.parametersSchema`
                        //           when this input is empty it will automatically use the configured schema or
                        //           "${projectDir}/nextflow_schema.json" as default. This input should not be empty
                        //           for meta pipelines
    help                // boolean:  show help message
    help_full           // boolean:  show full help message
    show_hidden         // boolean:  show hidden parameters in help message
    before_text         // string:   text to show before the help message and parameters summary
    after_text          // string:   text to show after the help message and parameters summary
    command             // string:   an example command of the pipeline

    main:

    if(help || help_full) {
        help_options = [
            beforeText: before_text,
            afterText: after_text,
            command: command,
            showHidden: show_hidden,
            fullHelp: help_full,
        ]
        if(parameters_schema) {
            help_options << [parametersSchema: parameters_schema]
        }
        log.info paramsHelp(
            help_options,
            params.help instanceof String ? params.help : "",
        )
        exit 0
    }

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //

    summary_options = [:]
    if(parameters_schema) {
        summary_options << [parametersSchema: parameters_schema]
    }
    log.info before_text
    log.info paramsSummaryLog(summary_options, input_workflow)
    log.info after_text

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(validate_params) {
        validateOptions = [:]
        if(parameters_schema) {
            validateOptions << [parametersSchema: parameters_schema]
        }
        validateParameters(validateOptions)
    }

    emit:
    dummy_emit = true
}

