import org.yaml.snakeyaml.Yaml

/*
 * This file holds several functions used to perform standard checks for the nf-core pipeline template.
 */

class Checks {

    static void check_conda_channels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // check that all channels are present
        def required_channels = ['conda-forge', 'bioconda', 'defaults']
        def conda_check_failed = !required_channels.every { ch -> ch in channels }

        // check that they are in the right order
        conda_check_failed |= !(channels.indexOf('conda-forge') < channels.indexOf('bioconda'))
        conda_check_failed |= !(channels.indexOf('bioconda') < channels.indexOf('defaults'))

        if (conda_check_failed) {
            log.warn(
                "There's a problem with your conda configuration. " +
                "Please refer to  https://bioconda.github.io/user/install.html#set-up-channels " + 
                "to learn how to set-up the conda-forge and bioconda channels correctly. " + 
                "NB: the order of the channels matters"
            )
        }
    }

    static void aws_batch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            assert (params.awsqueue && params.awsregion) : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            // related: https://github.com/nextflow-io/nextflow/issues/813
            assert params.outdir.startsWith('s3:')       : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
            // Prevent trace files to be stored on S3 since S3 does not support rolling files.
            assert !params.tracedir.startsWith('s3:')    :  "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
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

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.1400710\n\n" +
               "* The nf-core framework\n" +
               "  https://dx.doi.org/10.1038/s41587-020-0439-x\n" +
               "  https://rdcu.be/b1GjZ\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    // Exit pipeline if incorrect --genome key provided
    static void genome_exists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "============================================================================="
            System.exit(0)
        }
    }

    // Get attribute from genome config file e.g. fasta
    static String get_genome_attribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }

    // Print a warning after SRA download has completed
    static void sra_download(log) {
        log.warn "=============================================================================\n" +
                 "  THIS IS AN EXPERIMENTAL FEATURE!\n\n" + 
                 "  Please double-check the samplesheet that has been auto-created using the\n" +
                 "  public database ids provided via the '--public_data_ids' parameter.\n\n" +
                 "  Public databases don't reliably hold information such as experimental group,\n" +
                 "  replicate identifiers or strandedness information.\n\n" +  
                 "  All of the sample metadata obtained from the ENA has been appended\n" +
                 "  as additional columns to help you manually curate the samplesheet before\n" +
                 "  you run the pipeline.\n" +
                 "==================================================================================="
    }

    // Print a warning if using GRCh38 assembly from igenomes.config
    static void ncbi_genome_warn(log) {
        log.warn "=============================================================================\n" +
                 "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
                 "  Auto-activating '--skip_biotype_qc' parameter to circumvent the issue below:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/460.\n\n" +
                 "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
                 "==================================================================================="
    }

    // Print a warning if using a UCSC assembly from igenomes.config
    static void ucsc_genome_warn(log) {
        log.warn "=============================================================================\n" +
                 "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
                 "  Auto-activating '--skip_biotype_qc' parameter to circumvent the issue below:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/460.\n\n" +
                 "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
                 "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
                 "==================================================================================="
    }

    // Print a warning if both GTF and GFF have been provided
    static void gtf_gff_warn(log) {
        log.warn "=============================================================================\n" +
                 "  Both '--gtf' and '--gff' parameters have been provided.\n" +
                 "  Using GTF file as priority.\n" +
                 "==================================================================================="
    }

    // Print a warning if --skip_alignment has been provided
    static void skip_alignment_warn(log) {
        log.warn "=============================================================================\n" +
                 "  '--skip_alignment' parameter has been provided.\n" +
                 "  Skipping alignment, quantification and all downstream QC processes.\n" +
                 "==================================================================================="
    }

    // Print a warning if using '--aligner star_rsem' and '--with_umi'
    static void rsem_umi_error(log) {
        log.error "=============================================================================\n" +
                  "  When using '--aligner star_rsem', STAR is run by RSEM itself and so it is\n" +
                  "  not possible to remove UMIs before the quantification.\n\n" +
                  "  If you would like to remove UMI barcodes using the '--with_umi' option\n" + 
                  "  please use either '--aligner star_salmon' or '--aligner hisat2'.\n" +
                  "============================================================================="
        System.exit(0)
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
