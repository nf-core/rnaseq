/*
 * Functions to be run on completion of pipeline
 */

class Completion {
    static void email(workflow, params, summary_params, projectDir, log, multiqc_report=[], fail_percent_mapped=[:]) {

        // Set up the e-mail variables
        def subject = "[$workflow.manifest.name] Successful: $workflow.runName"
        if (fail_percent_mapped.size() > 0) {
            subject = "[$workflow.manifest.name] Partially successful (${fail_percent_mapped.size()} skipped): $workflow.runName"
        }
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
        email_fields['version']             = workflow.manifest.version
        email_fields['runName']             = workflow.runName
        email_fields['success']             = workflow.success
        email_fields['dateComplete']        = workflow.complete
        email_fields['duration']            = workflow.duration
        email_fields['exitStatus']          = workflow.exitStatus
        email_fields['errorMessage']        = (workflow.errorMessage ?: 'None')
        email_fields['errorReport']         = (workflow.errorReport ?: 'None')
        email_fields['commandLine']         = workflow.commandLine
        email_fields['projectDir']          = workflow.projectDir
        email_fields['summary']             = summary << misc_fields
        email_fields['fail_percent_mapped'] = fail_percent_mapped.keySet()
        email_fields['min_mapped_reads']    = params.min_mapped_reads
        
        // On success try attach the multiqc report
        def mqc_report = null
        try {
            if (workflow.success && !params.skip_multiqc) {
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
        def max_multiqc_email_size = params.max_multiqc_email_size as nextflow.util.MemoryUnit 
        def smail_fields           = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize:  max_multiqc_email_size.toBytes()]
        def sf                     = new File("$projectDir/assets/sendmail_template.txt")
        def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
        def sendmail_html          = sendmail_template.toString()

        // Send the HTML e-mail
        Map colors = Utils.logColours(params.monochrome_logs)
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

    static void summary(workflow, params, log, fail_percent_mapped=[:], pass_percent_mapped=[:]) {
        Map colors = Utils.logColours(params.monochrome_logs)

        def total_aln_count = pass_percent_mapped.size() + fail_percent_mapped.size()
        if (pass_percent_mapped.size() > 0) {
            def idx = 0
            def samp_aln = ''
            for (samp in pass_percent_mapped) {
                samp_aln += "    ${samp.value}%: ${samp.key}\n"
                idx += 1
                if (idx > 5) {
                    samp_aln += "    ..see pipeline reports for full list\n"
                    break;
                }
            }
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} ${pass_percent_mapped.size()}/$total_aln_count samples passed STAR ${params.min_mapped_reads}% mapped threshold:\n${samp_aln}${colors.reset}-"
        }
        if (fail_percent_mapped.size() > 0) {
            def samp_aln = ''
            for (samp in fail_percent_mapped) {
                samp_aln += "    ${samp.value}%: ${samp.key}\n"
            }
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} ${fail_percent_mapped.size()}/$total_aln_count samples skipped since they failed STAR ${params.min_mapped_reads}% mapped threshold:\n${samp_aln}${colors.reset}-"
        }

        if (workflow.success) {
            if (workflow.stats.ignoredCount == 0) {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
            } else {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
            }
        } else {
            Checks.hostName(workflow, params, log)
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
        }
    }
}
