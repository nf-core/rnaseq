/*
 * This file holds several functions used to perform standard checks for the nf-core pipeline template.
 */

class Checks {

    static void aws_batch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            assert !params.awsqueue || !params.awsregion : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            // related: https://github.com/nextflow-io/nextflow/issues/813
            assert !params.outdir.startsWith('s3:') : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
            // Prevent trace files to be stored on S3 since S3 does not support rolling files.
            assert params.tracedir.startsWith('s3:') :  "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
        }
    }

    static void hostname(workflow, params, log) {
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (params.hostnames) {
            def hostname = "hostname".execute().text.trim()
            params.hostnames.each { prof, hnames ->
                hnames.each { hname ->
                    if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                        log.info "=${colors.yellow}====================================================${colors.reset}=\n" +
                                  "${colors.yellow}WARN: You are running with `-profile $workflow.profile`\n" +
                                  "      but your machine hostname is ${colors.white}'$hostname'${colors.reset}.\n" +
                                  "      ${colors.yellow_bold}Please use `-profile $prof${colors.reset}`\n" +
                                  "=${colors.yellow}====================================================${colors.reset}="
                    }
                }
            }
        }
    }

    // Print a warning if using GRCh38 assembly from igenomes.config
    static void genome_warn(log) {
        log.warn "=============================================================================\n" +
                 "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
                 "  Auto-activating --skip_biotype_qc parameter to circumvent the issue below:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/460.\n\n" +
                 "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
                 "==================================================================================="
    }

    // Function that parses and returns the alignment rate from the STAR log output
    static Float get_star_percent_mapped(workflow, params, log, align_log) {
        def percent_aligned = 0
        def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        align_log.eachLine { line ->
            def matcher = line =~ pattern
            if (matcher) {
                percent_aligned = matcher[0][1].toFloat()
            }
        }
        def logname = align_log.getBaseName() - '.Log.final'
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (percent_aligned <= params.min_mapped_reads.toFloat()) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} [FAIL] STAR ${params.min_mapped_reads}% mapped threshold. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS: ${percent_aligned}% - $logname${colors.reset}."
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} [PASS] STAR ${params.min_mapped_reads}% mapped threshold: ${percent_aligned}% - $logname${colors.reset}."
        }
        return percent_aligned
    }

    static String fail_star_percent_mapped_multiqc(failed_samples) {
        String yaml_file_text  = """
        id: 'nf-core-rnaseq-fail-alignment'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/rnaseq Samples Failed Alignment'
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
