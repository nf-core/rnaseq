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

    // Print a warning if using '--aligner star_rsem' and '--with_umi'
    static void rsem_umi_error(log) {
        log.error "=============================================================================\n" +
                  "  When using '--aligner star_rsem', STAR is run by RSEM itself and so it is\n" +
                  "  not possible to remove UMIs before the quantification.\n\n" +
                  "  If you would like to remove UMI barcodes using the '--with_umi' option\n" + 
                  "  please use either '--aligner star' or '--aligner hisat2'.\n" +
                  "============================================================================="
    }

    // Function that parses and returns the alignment rate from the STAR log output
    static ArrayList get_star_percent_mapped(workflow, params, log, align_log) {
        def percent_aligned = 0
        def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        align_log.eachLine { line ->
            def matcher = line =~ pattern
            if (matcher) {
                percent_aligned = matcher[0][1].toFloat()
            }
        }

        def pass = false
        def logname = align_log.getBaseName() - '.Log.final'
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (percent_aligned <= params.min_mapped_reads.toFloat()) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} [FAIL] STAR ${params.min_mapped_reads}% mapped threshold. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS: ${percent_aligned}% - $logname${colors.reset}."
        } else {
            pass = true
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} [PASS] STAR ${params.min_mapped_reads}% mapped threshold: ${percent_aligned}% - $logname${colors.reset}."
        }
        return [ percent_aligned, pass ]
    }

    // Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
    static ArrayList get_inferexperiment_strandedness(inferexperiment_file, cutoff=30) {
        def sense        = 0
        def antisense    = 0
        def undetermined = 0
        inferexperiment_file.eachLine { line ->
            def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
            def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
            def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
            def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
            def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
            if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
            if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
            if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
            if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
            if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
        }
        def strandedness = 'unstranded'
        if (sense >= 100-cutoff) {
            strandedness = 'forward'
        } else if (antisense >= 100-cutoff) {
            strandedness = 'reverse'
        }
        return [ strandedness, sense, antisense, undetermined ]
    }
}
