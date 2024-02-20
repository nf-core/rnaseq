//
// Subworkflow with utility functions specific to the nf-core pipeline template
//

import org.yaml.snakeyaml.Yaml
import nextflow.extension.FilesEx

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFCORE_PIPELINE {

    take:
    nextflow_cli_args

    main:
    valid_config = checkConfigProvided()
    checkProfileProvided(nextflow_cli_args)

    emit:
    valid_config
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
//  Warn if a -profile or Nextflow config has not been provided to run the pipeline
//
def checkConfigProvided() {
    valid_config = true
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
            "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
            "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
            "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
            "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
            "Please refer to the quick start section and usage docs for the pipeline.\n "
        valid_config = false
    }
    return valid_config
}

//
// Exit pipeline if --profile contains spaces
//
def checkProfileProvided(nextflow_cli_args) {
    if (workflow.profile.endsWith(',')) {
        error "The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
    if (nextflow_cli_args[0]) {
        log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${nextflow_cli_args[0]}` has been detected.\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
}

//
// Citation string for pipeline
//
def workflowCitation() {
    return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The pipeline\n" +
        "  ${workflow.manifest.doi}\n\n" +
        "* The nf-core framework\n" +
        "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}

//
// Generate workflow version string
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
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { processVersionsFromYAML(it) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}

//
// Get workflow summary for MultiQC
//
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    for (group in summary_params.keySet()) {
        def group_params = summary_params.get(group)  // This gets the parameters of that particular group
        if (group_params) {
            summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
            summary_section += "    <dl class=\"dl-horizontal\">\n"
            for (param in group_params.keySet()) {
                summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
            }
            summary_section += "    </dl>\n"
        }
    }

    String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
    yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text        += "plot_type: 'html'\n"
    yaml_file_text        += "data: |\n"
    yaml_file_text        += "${summary_section}"

    return yaml_file_text
}

//
// nf-core logo
//
def nfCoreLogo(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
                                                ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
        ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
        ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
        ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                ${colors.green}`._,._,\'${colors.reset}
        ${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}

//
// Return dashed line
//
def dashedLine(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

//
// ANSII colours used for terminal logging
//
def logColours(monochrome_logs=true) {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

//
// Attach the multiqc report to email
//
def attachMultiqcReport(multiqc_report) {
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList && mqc_report.size() >= 1) {
                if (mqc_report.size() > 1) {
                    log.warn "[$workflow.manifest.name] Found multiple reports from process 'MULTIQC', will use only one"
                }
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        if (multiqc_report) {
            log.warn "[$workflow.manifest.name] Could not attach MultiQC report to summary email"
        }
    }
    return mqc_report
}

//
// Construct and send completion email
//
def completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs=true, multiqc_report=null) {

    // Set up the e-mail variables
    def subject = "[$workflow.manifest.name] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[$workflow.manifest.name] FAILED: $workflow.runName"
    }

    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    def misc_fields = [:]
    misc_fields['Date Started']              = workflow.start
    misc_fields['Date Completed']            = workflow.complete
    misc_fields['Pipeline script file path'] = workflow.scriptFile
    misc_fields['Pipeline script hash ID']   = workflow.scriptId
    if (workflow.repository) misc_fields['Pipeline repository Git URL']    = workflow.repository
    if (workflow.commitId)   misc_fields['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision)   misc_fields['Pipeline Git branch/tag']        = workflow.revision
    misc_fields['Nextflow Version']           = workflow.nextflow.version
    misc_fields['Nextflow Build']             = workflow.nextflow.build
    misc_fields['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    def email_fields = [:]
    email_fields['version']      = getWorkflowVersion()
    email_fields['runName']      = workflow.runName
    email_fields['success']      = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration']     = workflow.duration
    email_fields['exitStatus']   = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport']  = (workflow.errorReport ?: 'None')
    email_fields['commandLine']  = workflow.commandLine
    email_fields['projectDir']   = workflow.projectDir
    email_fields['summary']      = summary << misc_fields

    // On success try attach the multiqc report
    def mqc_report = attachMultiqcReport(multiqc_report)

    // Check if we are only sending emails on failure
    def email_address = email
    if (!email && email_on_fail && !workflow.success) {
        email_address = email_on_fail
    }

    // Render the TXT template
    def engine       = new groovy.text.GStringTemplateEngine()
    def tf           = new File("${workflow.projectDir}/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt    = txt_template.toString()

    // Render the HTML template
    def hf            = new File("${workflow.projectDir}/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html    = html_template.toString()

    // Render the sendmail template
    def max_multiqc_email_size = (params.containsKey('max_multiqc_email_size') ? params.max_multiqc_email_size : 0) as nextflow.util.MemoryUnit
    def smail_fields           = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "${workflow.projectDir}", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes() ]
    def sf                     = new File("${workflow.projectDir}/assets/sendmail_template.txt")
    def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html          = sendmail_template.toString()

    // Send the HTML e-mail
    Map colors = logColours(monochrome_logs)
    if (email_address) {
        try {
            if (plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            def sendmail_tf = new File(workflow.launchDir.toString(), ".sendmail_tmp.html")
            sendmail_tf.withWriter { w -> w << sendmail_html }
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (sendmail)-"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (mail)-"
        }
    }

    // Write summary e-mail HTML to a file
    def output_hf = new File(workflow.launchDir.toString(), ".pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    FilesEx.copyTo(output_hf.toPath(), "${outdir}/pipeline_info/pipeline_report.html");
    output_hf.delete()

    // Write summary e-mail TXT to a file
    def output_tf = new File(workflow.launchDir.toString(), ".pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }
    FilesEx.copyTo(output_tf.toPath(), "${outdir}/pipeline_info/pipeline_report.txt");
    output_tf.delete()
}

//
// Print pipeline summary on completion
//
def completionSummary(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
        }
    } else {
        log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
    }
}

//
// Construct and send a notification to a web server as JSON e.g. Microsoft Teams and Slack
//
def imNotification(summary_params, hook_url) {
    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    def misc_fields = [:]
    misc_fields['start']                                = workflow.start
    misc_fields['complete']                             = workflow.complete
    misc_fields['scriptfile']                           = workflow.scriptFile
    misc_fields['scriptid']                             = workflow.scriptId
    if (workflow.repository) misc_fields['repository']  = workflow.repository
    if (workflow.commitId)   misc_fields['commitid']    = workflow.commitId
    if (workflow.revision)   misc_fields['revision']    = workflow.revision
    misc_fields['nxf_version']                          = workflow.nextflow.version
    misc_fields['nxf_build']                            = workflow.nextflow.build
    misc_fields['nxf_timestamp']                        = workflow.nextflow.timestamp

    def msg_fields = [:]
    msg_fields['version']      = getWorkflowVersion()
    msg_fields['runName']      = workflow.runName
    msg_fields['success']      = workflow.success
    msg_fields['dateComplete'] = workflow.complete
    msg_fields['duration']     = workflow.duration
    msg_fields['exitStatus']   = workflow.exitStatus
    msg_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    msg_fields['errorReport']  = (workflow.errorReport ?: 'None')
    msg_fields['commandLine']  = workflow.commandLine.replaceFirst(/ +--hook_url +[^ ]+/, "")
    msg_fields['projectDir']   = workflow.projectDir
    msg_fields['summary']      = summary << misc_fields

    // Render the JSON template
    def engine       = new groovy.text.GStringTemplateEngine()
    // Different JSON depending on the service provider
    // Defaults to "Adaptive Cards" (https://adaptivecards.io), except Slack which has its own format
    def json_path     = hook_url.contains("hooks.slack.com") ? "slackreport.json" : "adaptivecard.json"
    def hf            = new File("${workflow.projectDir}/assets/${json_path}")
    def json_template = engine.createTemplate(hf).make(msg_fields)
    def json_message  = json_template.toString()

    // POST
    def post = new URL(hook_url).openConnection();
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")
    post.getOutputStream().write(json_message.getBytes("UTF-8"));
    def postRC = post.getResponseCode();
    if (! postRC.equals(200)) {
        log.warn(post.getErrorStream().getText());
    }
}
