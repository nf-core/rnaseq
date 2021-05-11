/*
 * This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
 */

class WorkflowRnaseq {

    /*
     * Check and validate parameters
     */
    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.gtf && !params.gff) {
            log.error "No GTF or GFF3 annotation specified! The pipeline requires at least one of these files."
            System.exit(1)
        }

        if (params.gtf) {
            if (params.gff) {
                gtfGffWarn(log)
            }
            if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
                ncbiGenomeWarn(log)
            }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }

        if (!params.skip_alignment) {
            if (!valid_params['aligners'].contains(params.aligner)) {
                log.error "Invalid option: '${params.aligner}'. Valid options for '--aligner': ${valid_params['aligners'].join(', ')}."
                System.exit(1)
            }
        } else {
            if (!params.pseudo_aligner) {
                log.error "--skip_alignment specified without --pseudo_aligner...please specify e.g. --pseudo_aligner ${valid_params['pseudoaligners'][0]}."
                System.exit(1)
            }
            skipAlignmentWarn(log)
        }

        if (params.pseudo_aligner) {
            if (!valid_params['pseudoaligners'].contains(params.pseudo_aligner)) {
                log.error "Invalid option: '${params.pseudo_aligner}'. Valid options for '--pseudo_aligner': ${valid_params['pseudoaligners'].join(', ')}."
                System.exit(1)
            } else {
                if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
                    log.error "To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff."
                    System.exit(1)
                }
            }
        }

        // Checks when running --aligner star_rsem
        if (!params.skip_alignment && params.aligner == 'star_rsem') {
            if (params.with_umi) {
                rsemUmiError(log)
            }
            if (params.rsem_index && params.star_index) {
                rsemStarIndexWarn(log)
            }
        }

        // Check which RSeQC modules we are running
        def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['rseqc_modules'] + rseqc_modules).unique().size() != valid_params['rseqc_modules'].size()) {
            log.error "Invalid option: ${params.rseqc_modules}. Valid options for '--rseqc_modules': ${valid_params['rseqc_modules'].join(', ')}"
            System.exit(1)
        }
    }

    /*
     * Function to check whether biotype field exists in GTF file
     */
    public static Boolean biotypeInGtf(gtf_file, biotype, log) {
        def hits = 0
        gtf_file.eachLine { line ->
            def attributes = line.split('\t')[-1].split()
            if (attributes.contains(biotype)) {
                hits += 1
            }
        }
        if (hits) {
            return true
        } else {
            log.warn "=============================================================================\n" +
                "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
                "  Biotype QC will be skipped to circumvent the issue below:\n" +
                "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
                "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
                "==================================================================================="
            return false
        }
    }

    /*
     * Function that parses and returns the alignment rate from the STAR log output
     */
    public static ArrayList getStarPercentMapped(params, align_log) {
        def percent_aligned = 0
        def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        align_log.eachLine { line ->
            def matcher = line =~ pattern
            if (matcher) {
                percent_aligned = matcher[0][1].toFloat()
            }
        }

        def pass = false
        if (percent_aligned >= params.min_mapped_reads.toFloat()) {
            pass = true
        }
        return [ percent_aligned, pass ]
    }

    /*
     * Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
     */
    public static ArrayList getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
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

    /*
     * Get workflow summary for MultiQC
     */
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
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

    /*
     * Exit pipeline if incorrect --genome key provided
     */
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }

    /*
     * Print a warning if using GRCh38 assembly from igenomes.config
     */
    private static void ncbiGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    /*
     * Print a warning if using a UCSC assembly from igenomes.config
     */
    private static void ucscGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    /*
     * Print a warning if both GTF and GFF have been provided
     */
    private static void gtfGffWarn(log) {
        log.warn "=============================================================================\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "==================================================================================="
    }

    /*
     * Print a warning if --skip_alignment has been provided
     */
    private static void skipAlignmentWarn(log) {
        log.warn "=============================================================================\n" +
            "  '--skip_alignment' parameter has been provided.\n" +
            "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
            "==================================================================================="
    }

    /*
     * Print a warning if using '--aligner star_rsem' and '--with_umi'
     */
    private static void rsemUmiError(log) {
        log.error "=============================================================================\n" +
            "  When using '--aligner star_rsem', STAR is run by RSEM itself and so it is\n" +
            "  not possible to remove UMIs before the quantification.\n\n" +
            "  If you would like to remove UMI barcodes using the '--with_umi' option\n" +
            "  please use either '--aligner star_salmon' or '--aligner hisat2'.\n" +
            "============================================================================="
        System.exit(1)
    }

    /*
     * Print a warning if using '--aligner star_rsem' and providing both '--rsem_index' and '--star_index'
     */
    private static void rsemStarIndexWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using '--aligner star_rsem', both the STAR and RSEM indices should\n" +
            "  be present in the path specified by '--rsem_index'.\n\n" +
            "  This warning has been generated because you have provided both\n" +
            "  '--rsem_index' and '--star_index'. The pipeline will ignore the latter.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/568\n" +
            "==================================================================================="
    }
}
