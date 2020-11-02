/*
 * This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
 */

import groovy.json.JsonSlurper

class Schema {
    /*
     * This method tries to read a JSON params file
     */
    private static LinkedHashMap params_load(String schema_path) {
        def params_map = new LinkedHashMap()
        try {
            params_map = params_read(schema_path)
        } catch (Exception e) {
            println "Could not read parameters settings from JSON. $e"
            params_map = new LinkedHashMap()
        }
        return params_map
    }

    /*
    Method to actually read in JSON file using Groovy.
    Group (as Key), values are all parameters
        - Parameter1 as Key, Description as Value
        - Parameter2 as Key, Description as Value
        ....
    Group
        -
    */
    private static LinkedHashMap params_read(String schema_path) throws Exception {
        def json = new File(schema_path).text
        def Map json_params = (Map) new JsonSlurper().parseText(json).get('definitions')
        /* Tree looks like this in nf-core schema
         * definitions <- this is what the first get('definitions') gets us
             group 1
               title
               description
                 properties
                   parameter 1
                     type
                     description
                   parameter 2
                     type
                     description
             group 2
               title
               description
                 properties
                   parameter 1
                     type
                     description
        */
        def params_map = new LinkedHashMap()
        json_params.each { key, val ->
            def Map group = json_params."$key".properties // Gets the property object of the group
            def title = json_params."$key".title
            def sub_params = new LinkedHashMap()
            group.each { innerkey, value ->
                sub_params.put(innerkey, value)
            }
            params_map.put(title, sub_params)
        }
        return params_map
    }

    /*
     * Get maximum number of characters across all parameter names
     */
    private static Integer params_max_chars(params_map) {
        Integer max_chars = 0
        for (group in params_map.keySet()) {
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            for (param in group_params.keySet()) {
                if (param.size() > max_chars) {
                    max_chars = param.size()
                }
            }
        }
        return max_chars
    }

    /*
     * Beautify parameters for --help
     */
    private static String params_help_string(schema_path, command, monochrome_logs) {
        String output = "Typical pipeline command:\n\n"
        output += "    ${command}\n\n"
        def params_map = params_load(schema_path)
        def max_chars = params_max_chars(params_map) + 1
        for (group in params_map.keySet()) {
            output += group + "\n"
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            for (param in group_params.keySet()) {
                def type = "[" + group_params.get(param).type + "]"
                def description = group_params.get(param).description
                output+= "    \u001B[1m" +  param.padRight(max_chars) + "\u001B[1m" + type.padRight(10) + description + "\n"
            }
            output += "\n"
        }
        output += Headers.dashed_line(monochrome_logs)
        return output
    }

    /*
     * Groovy Map summarising parameters/workflow options used by the pipeline
     */
    private static LinkedHashMap params_summary_map(workflow, params, schema_path, force_params=[]) {
        // Get pipeline parameters defined in JSON Schema
        def Map params_summary = [:]
        def params_map = params_load(schema_path)
        for (group in params_map.keySet()) {
            def sub_params = new LinkedHashMap()
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            for (param in group_params.keySet()) {
                if (params.containsKey(param)) {
                    def value = params.get(param)
                    if (!group_params.get(param).hidden || force_params.contains(param)) {
                        sub_params.put("$param", value)
                    }
                }
            }
            params_summary.put(group, sub_params)
        }

        // Get a selection of core Nextflow workflow options
        def Map workflow_summary = [:]        
        if (workflow.revision) {
            workflow_summary['revision'] = workflow.revision
        }
        workflow_summary['runName']      = workflow.runName
        if (workflow.containerEngine) {
            workflow_summary['containerEngine'] = "$workflow.containerEngine"
            workflow_summary['container']       = "$workflow.container"
        }
        workflow_summary['launchDir']    = workflow.launchDir
        workflow_summary['workDir']      = workflow.workDir
        workflow_summary['projectDir']   = workflow.projectDir
        workflow_summary['userName']     = workflow.userName
        workflow_summary['profile']      = workflow.profile
        workflow_summary['configFiles']  = workflow.configFiles.join(', ')
        
        return [ 'Nextflow workflow options' : workflow_summary ] << params_summary
    }

    /*
     * Beautify parameters for summary and return as string
     */
    private static String params_summary_string(params_map, monochrome_logs) {
        String output = ""
        def max_chars = params_max_chars(params_map)
        for (group in params_map.keySet()) {
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                output += group + "\n"
                for (param in group_params.keySet()) {
                    output+= "    \u001B[1m" +  param.padRight(max_chars) + ": \u001B[1m" + group_params.get(param) + "\n"
                }
                output += "\n"
            }
        }
        output += Headers.dashed_line(monochrome_logs)
        return output
    }

    static String params_summary_multiqc(summary) {
        String yaml_file_text  = """
        id: 'nf-core-rnaseq-summary'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/rnaseq Workflow Summary'
        section_href: 'https://github.com/nf-core/rnaseq'
        plot_type: 'html'
        data: |
            <dl class=\"dl-horizontal\">
            ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
            </dl>
        """.stripIndent()
        return yaml_file_text
    }
}
