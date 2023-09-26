//
// This file holds several functions used within the nf-core pipeline template.
//

import org.yaml.snakeyaml.Yaml
import groovy.json.JsonOutput

class NfcoreTemplate {

    //
    // Check AWS Batch related parameters have been specified correctly
    //
    public static void awsBatch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            // Check params.awsqueue and params.awsregion have been set if running on AWSBatch
            assert (params.awsqueue && params.awsregion) : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            assert params.outdir.startsWith('s3:')       : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
        }
    }

    //
    //  Warn if a -profile or Nextflow config has not been provided to run the pipeline
    //
    public static void checkConfigProvided(workflow, log) {
        if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
            log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
                    "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
                    "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
                    "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
                    "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
                    "Please refer to the quick start section and usage docs for the pipeline.\n "
        }
    }

    //
    // Generate version string
    //
    public static String version(workflow) {
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
    // Construct and send completion email
    //
    public static void email(workflow, params, summary_params, projectDir, log, multiqc_report=[]) {

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
        email_fields['version']      = NfcoreTemplate.version(workflow)
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

        // Check if we are only sending emails on failure
        def email_address = params.email
        if (!params.email && params.email_on_fail && !workflow.success) {
            email_address = params.email_on_fail
        }

        // Render the TXT template
        def engine       = new groovy.text.GStringTemplateEngine()
        def tf           = new File("$projectDir/assets/email_template.txt")
        def txt_template = engine.createTemplate(tf).make(email_fields)
        def email_txt    = txt_template.toString()

        // Render the HTML template
        def hf            = new File("$projectDir/assets/email_template.html")
        def html_template = engine.createTemplate(hf).make(email_fields)
        def email_html    = html_template.toString()

        // Render the sendmail template
        def max_multiqc_email_size = (params.containsKey('max_multiqc_email_size') ? params.max_multiqc_email_size : 0) as nextflow.util.MemoryUnit
        def smail_fields           = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes() ]
        def sf                     = new File("$projectDir/assets/sendmail_template.txt")
        def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
        def sendmail_html          = sendmail_template.toString()

        // Send the HTML e-mail
        Map colors = logColours(params.monochrome_logs)
        if (email_address) {
            try {
                if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
                // Try to send HTML e-mail using sendmail
                [ 'sendmail', '-t' ].execute() << sendmail_html
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (sendmail)-"
            } catch (all) {
                // Catch failures and try with plaintext
                def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
                if ( mqc_report.size() <= max_multiqc_email_size.toBytes() ) {
                    mail_cmd += [ '-A', mqc_report ]
                }
                mail_cmd.execute() << email_html
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (mail)-"
            }
        }

        // Write summary e-mail HTML to a file
        def output_d = new File("${params.outdir}/pipeline_info/")
        if (!output_d.exists()) {
            output_d.mkdirs()
        }
        def output_hf = new File(output_d, "pipeline_report.html")
        output_hf.withWriter { w -> w << email_html }
        def output_tf = new File(output_d, "pipeline_report.txt")
        output_tf.withWriter { w -> w << email_txt }
    }

    //
    // Construct and send a notification to a web server as JSON
    // e.g. Microsoft Teams and Slack
    //
    public static void IM_notification(workflow, params, summary_params, projectDir, log) {
        def hook_url = params.hook_url

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
        msg_fields['version']      = NfcoreTemplate.version(workflow)
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
        def hf            = new File("$projectDir/assets/${json_path}")
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

    //
    // Dump pipeline parameters in a json file
    //
    public static void dump_parameters(workflow, params) {
        def output_d = new File("${params.outdir}/pipeline_info/")
        if (!output_d.exists()) {
            output_d.mkdirs()
        }

        def timestamp  = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
        def output_pf  = new File(output_d, "params_${timestamp}.json")
        def jsonStr    = JsonOutput.toJson(params)
        output_pf.text = JsonOutput.prettyPrint(jsonStr)
    }

    //
    // Print pipeline summary on completion
    //
    public static void summary(workflow, params, log) {
        Map colors = logColours(params.monochrome_logs)
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
    // ANSII Colours used for terminal logging
    //
    public static Map logColours(Boolean monochrome_logs) {
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
    // Does what is says on the tin
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }

    //
    // nf-core logo
    //
    public static String logo(workflow, monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        String workflow_version = NfcoreTemplate.version(workflow)
        String.format(
            """\n
            ${dashedLine(monochrome_logs)}
                                                    ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
            ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
            ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
            ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                    ${colors.green}`._,._,\'${colors.reset}
            ${colors.purple}  ${workflow.manifest.name} ${workflow_version}${colors.reset}
            ${dashedLine(monochrome_logs)}
            """.stripIndent()
        )
    }
}
