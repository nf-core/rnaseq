//
// Subworkflow with functionality that may be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml
import groovy.json.JsonOutput
import nextflow.extension.FilesEx

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NEXTFLOW_PIPELINE {

    take:
    print_version        // boolean: print version
    dump_parameters      // boolean: dump parameters
    outdir               //    path: base directory used to publish pipeline results
    check_conda_channels // boolean: check conda channels

    main:

    //
    // Print workflow version and exit on --version
    //
    if (print_version) {
        log.info "${workflow.manifest.name} ${getWorkflowVersion()}"
        System.exit(0)
    }

    //
    // Dump pipeline parameters to a JSON file
    //
    if (dump_parameters && outdir) {
        dumpParametersToJSON(outdir)
    }

    //
    // When running with Conda, warn if channels have not been set-up appropriately
    //
    if (check_conda_channels) {
        checkCondaChannels()
    }

    emit:
    dummy_emit = true
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Generate version string
//
def getWorkflowVersion() {
    String version_string = ""
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Dump pipeline parameters to a JSON file
//
def dumpParametersToJSON(outdir) {
    def timestamp  = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
    def filename   = "params_${timestamp}.json"
    def temp_pf    = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr    = JsonOutput.toJson(params)
    temp_pf.text   = JsonOutput.prettyPrint(jsonStr)

    FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
}

//
// When running with -profile conda, warn if channels have not been set-up appropriately
//
def checkCondaChannels() {
    Yaml parser = new Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    } catch(NullPointerException | IOException e) {
        log.warn "Could not verify conda channel configuration."
        return
    }

    // Check that all channels are present
    // This channel list is ordered by required channel priority.
    def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
    def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

    // Check that they are in the right order
    def channel_priority_violation = false
    def n = required_channels_in_order.size()
    for (int i = 0; i < n - 1; i++) {
        channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
    }

    if (channels_missing | channel_priority_violation) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  There is a problem with your Conda configuration!\n\n" +
            "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
            "  Please refer to https://bioconda.github.io/\n" +
            "  The observed channel order is \n" +
            "  ${channels}\n" +
            "  but the following channel order is required:\n" +
            "  ${required_channels_in_order}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }
}
