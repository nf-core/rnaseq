//
// Subworkflow with functionality specific to the nf-core/rnaseq pipeline
//

import groovy.json.JsonSlurper

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { calculateStrandedness     } from '../../nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    emit:
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report
    trim_status        // map: pass/fail status per sample for trimming
    map_status         // map: pass/fail status per sample for mapping
    strand_status      // map: pass/fail status per sample for strandedness check

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    trim_status
        .map{
            id, status -> pass_trimmed_reads[id] = status
        }

    map_status
        .map{
            id, status -> pass_mapped_reads[id] = status
        }

    strand_status
        .map{
            id, status -> pass_strand_check[id] = status
        }

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList()
            )
        }

        rnaseqSummary(monochrome_logs=monochrome_logs, pass_mapped_reads=pass_mapped_reads, pass_trimmed_reads=pass_trimmed_reads, pass_strand_check=pass_strand_check)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Function to check samples are internally consistent after being grouped
//
def checkSamplesAfterGrouping(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect{ it.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    genomeExistsError()

    if (!params.fasta) {
        error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
    }

    if (!params.gtf && !params.gff) {
        error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
    }

    if (params.gtf) {
        if (params.gff) {
            gtfGffWarn()
        }
        if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
            ncbiGenomeWarn()
        }
        if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
            ucscGenomeWarn()
        }
    }

    if (params.transcript_fasta) {
        transcriptsFastaWarn()
    }

    if (!params.skip_bbsplit && !params.bbsplit_index && !params.bbsplit_fasta_list) {
        error("Please provide either --bbsplit_fasta_list / --bbsplit_index to run BBSplit.")
    }

    if (params.remove_ribo_rna && !params.ribo_database_manifest) {
        error("Please provide --ribo_database_manifest to remove ribosomal RNA with SortMeRNA.")
    }

    if (params.with_umi && !params.skip_umi_extract) {
        if (!params.umitools_bc_pattern && !params.umitools_bc_pattern2) {
            error("UMI-tools requires a barcode pattern to extract barcodes from the reads.")
        }
    }

    if (params.skip_alignment) {
        skipAlignmentWarn()
    }

    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {
        if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
            error("To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff.")
        }
    }

    // Checks when running --aligner star_rsem
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        if (params.with_umi) {
            rsemUmiError()
        }
        if (params.rsem_index && params.star_index) {
            rsemStarIndexWarn()
        }
        if (params.aligner  == 'star_rsem' && params.extra_star_align_args) {
            rsemStarExtraArgumentsWarn()
        }
    }

    // Warn if --additional_fasta provided with aligner index
    if (!params.skip_alignment && params.additional_fasta) {
        def index = ''
        if (params.aligner == 'star_salmon' && params.star_index) {
            index = 'star'
        }
        if (params.aligner == 'star_rsem' && params.rsem_index) {
            index = 'rsem'
        }
        if (params.aligner == 'hisat2' && params.hisat2_index) {
            index = 'hisat2'
        }
        if (index) {
            additionaFastaIndexWarn(index)
        }
    }

    //General checks for if contaminant screening is used
    if (params.contaminant_screening) {
        if (params.aligner == 'star_rsem') {
            error("Contaminant screening cannot be done with --aligner star_rsem since unaligned reads are not saved. Please use --aligner star_salmon or --aligner hisat2.")
        }
    }

    // Check that Kraken/Bracken database provided if using kraken2/bracken
    if (params.contaminant_screening in ['kraken2', 'kraken2_bracken']) {
        if (!params.kraken_db) {
            error("Contaminant screening set to kraken2 but not database is provided. Please provide a database with the --kraken_db option.")
        }
    // Check that Kraken/Bracken parameters are not provided when Kraken2 is not being used
    } else {
        if (!params.bracken_precision.equals('S')) {
            brackenPrecisionWithoutKrakenDBWarn()
        }

        if (params.save_kraken_assignments || params.save_kraken_unassigned || params.kraken_db) {
            krakenArgumentsWithoutKrakenDBWarn()
        }
    }

    // Check which RSeQC modules we are running
    def valid_rseqc_modules = ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
    def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
    if ((valid_rseqc_modules + rseqc_modules).unique().size() != valid_rseqc_modules.size()) {
        error("Invalid option: ${params.rseqc_modules}. Valid options for '--rseqc_modules': ${valid_rseqc_modules.join(', ')}")
    }

    // Check rRNA databases for sortmerna
    if (params.remove_ribo_rna) {
        ch_ribo_db = file(params.ribo_database_manifest)
        if (ch_ribo_db.isEmpty()) {
            error("File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!")
        }
    }

    // Check if file with list of fastas is provided when running BBSplit
    if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
        ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list)
        if (ch_bbsplit_fasta_list.isEmpty()) {
            error("File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!")
        }
    }
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()
    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Print a warning if both GTF and GFF have been provided
//
def gtfGffWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  Both '--gtf' and '--gff' parameters have been provided.\n" +
        "  Using GTF file as priority.\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using GRCh38 assembly from igenomes.config
//
def ncbiGenomeWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
        "  Biotype QC will be skipped to circumvent the issue below:\n" +
        "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
        "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using a UCSC assembly from igenomes.config
//
def ucscGenomeWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
        "  Biotype QC will be skipped to circumvent the issue below:\n" +
        "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
        "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using '--transcript_fasta'
//
def transcriptsFastaWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  '--transcript_fasta' parameter has been provided.\n" +
        "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
        "  Please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/753\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if --skip_alignment has been provided
//
def skipAlignmentWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  '--skip_alignment' parameter has been provided.\n" +
        "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using '--aligner star_rsem' and '--with_umi'
//
def rsemUmiError() {
    def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  When using '--aligner star_rsem', STAR is run by RSEM itself and so it is\n" +
        "  not possible to remove UMIs before the quantification.\n\n" +
        "  If you would like to remove UMI barcodes using the '--with_umi' option\n" +
        "  please use either '--aligner star_salmon' or '--aligner hisat2'.\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    error(error_string)
}

//
// Print a warning if using '--aligner star_rsem' and providing both '--rsem_index' and '--star_index'
//
def rsemStarIndexWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  When using '--aligner star_rsem', both the STAR and RSEM indices should\n" +
        "  be present in the path specified by '--rsem_index'.\n\n" +
        "  This warning has been generated because you have provided both\n" +
        "  '--rsem_index' and '--star_index'. The pipeline will ignore the latter.\n\n" +
        "  Please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/568\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using '--aligner star_rsem' and providing '--star_extra_alignment_args'
//
def rsemStarExtraArgumentsWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  No additional arguments can be passed to STAR when using RSEM.\n" +
        "  Because RSEM enforces its own parameters for STAR, any extra arguments\n" +
        "  to STAR will be ignored. Alternatively, choose the STAR+Salmon route.\n\n" +
        "  This warning has been generated because you have provided both\n" +
        "  '--aligner star_rsem' and '--extra_star_align_args'.\n\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if using '--additional_fasta' and '--<ALIGNER>_index'
//
def additionaFastaIndexWarn(index) {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  When using '--additional_fasta <FASTA_FILE>' the aligner index will not\n" +
        "  be re-built with the transgenes incorporated by default since you have \n" +
        "  already provided an index via '--${index}_index <INDEX>'.\n\n" +
        "  Set '--additional_fasta <FASTA_FILE> --${index}_index false --gene_bed false --save_reference'\n" +
        "  to re-build the index with transgenes included and the index and gene BED file will be saved in\n" +
        "  'results/genome/index/${index}/' for re-use with '--${index}_index'.\n\n" +
        "  Ignore this warning if you know that the index already contains transgenes.\n\n" +
        "  Please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/556\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Print a warning if --save_kraken_assignments or --save_kraken_unassigned is provided without --kraken_db
//
def krakenArgumentsWithoutKrakenDBWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  'Kraken2 related arguments have been provided without setting contaminant\n" +
        "  screening to Kraken2. Kraken2 is not being run so these will not be used.\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

///
/// Print a warning if --bracken-precision is provided without --kraken_db
///
def brackenPrecisionWithoutKrakenDBWarn() {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  '--bracken-precision' parameter has been provided without Kraken2 contaminant screening.\n" +
        "  Bracken will not run so precision will not be set.\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

//
// Function to generate an error if contigs in genome fasta file > 512 Mbp
//
def checkMaxContigSize(fai_file) {
    def max_size = 512000000
    fai_file.eachLine { line ->
        def lspl  = line.split('\t')
        def chrom = lspl[0]
        def size  = lspl[1]
        if (size.toInteger() > max_size) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
                "  ${chrom}: ${size}\n\n" +
                "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
                "  Please see:\n" +
                "  https://github.com/nf-core/rnaseq/issues/744\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            error(error_string)
        }
    }
}

//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(params, align_log) {
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

//
// Function to check whether biotype field exists in GTF file
//
def biotypeInGtf(gtf_file, biotype) {
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
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return false
    }
}

//
// Function that parses RSeQC infer_experiment output file to get inferred strandedness
//
def getInferexperimentStrandedness(inferexperiment_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def forwardFragments = 0
    def reverseFragments = 0
    def unstrandedFragments = 0

    inferexperiment_file.eachLine { line ->
        def unstranded_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
        def se_sense_matcher = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
        def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
        def pe_sense_matcher = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
        def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/

        if (unstranded_matcher) unstrandedFragments = unstranded_matcher[0][1].toFloat() * 100
        if (se_sense_matcher) forwardFragments = se_sense_matcher[0][1].toFloat() * 100
        if (se_antisense_matcher) reverseFragments = se_antisense_matcher[0][1].toFloat() * 100
        if (pe_sense_matcher) forwardFragments = pe_sense_matcher[0][1].toFloat() * 100
        if (pe_antisense_matcher) reverseFragments = pe_antisense_matcher[0][1].toFloat() * 100
    }

    // Use shared calculation function to determine strandedness
    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

//
// Print pipeline summary on completion
//
def rnaseqSummary(monochrome_logs=true, pass_mapped_reads=[:], pass_trimmed_reads=[:], pass_strand_check=[:]) {
    def colors = logColours(monochrome_logs)

    def fail_mapped_count  = pass_mapped_reads.count  { key, value -> value == false }
    def fail_trimmed_count = pass_trimmed_reads.count { key, value -> value == false }
    def fail_strand_count  = pass_strand_check.count  { key, value -> value == false }
    if (workflow.success) {
        def color = colors.green
        def status = []
        if (workflow.stats.ignoredCount != 0) {
            color = colors.yellow
            status += ['with errored process(es)']
        }
        if (fail_mapped_count > 0 || fail_trimmed_count > 0) {
            color = colors.yellow
            status += ['with skipped sampl(es)']
        }
        log.info "-${colors.purple}[$workflow.manifest.name]${color} Pipeline completed successfully ${status.join(', ')}${colors.reset}-"
        if (fail_trimmed_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${fail_trimmed_count}/${pass_trimmed_reads.size()} samples skipped since they failed ${params.min_trimmed_reads} trimmed read threshold.${colors.reset}-"
        }
        if (fail_mapped_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${fail_mapped_count}/${pass_mapped_reads.size()} samples skipped since they failed STAR ${params.min_mapped_reads}% mapped threshold.${colors.reset}-"
        }
        if (fail_strand_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${fail_strand_count}/${pass_strand_check.size()} samples failed strandedness check.${colors.reset}-"
        }
    } else {
        log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
    }
}

