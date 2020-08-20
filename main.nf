#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnaseq
========================================================================================
 nf-core/rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnaseq --input '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker

    Mandatory arguments:
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Generic:
      --single_end [bool]             Specifies that the input is single-end reads
      --sample_level [bool]           Used to turn off the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples

    References:                       If not specified in the configuration file or you wish to overwrite any of the references.
      --genome [str]                  Name of iGenomes reference
      --star_index [file]             Path to STAR index
      --hisat2_index [file]           Path to HiSAT2 index
      --salmon_index [file]           Path to Salmon index
      --fasta [file]                  Path to genome fasta file
      --transcript_fasta [file]       Path to transcript fasta file
      --additional_fasta [file]       Additional fasta file(s) containing e.g. ERCCs spike-ins, GFP or CAR-T transgene sequences to map to
      --splicesites [file]            Path to splice sites file for building HiSat2 index
      --gtf [file]                    Path to GTF file
      --gff [file]                    Path to GFF3 file
      --bed12 [file]                  Path to bed12 file
      --star_index_options [str]      Additional options that will be appended to the STAR genome indexing command
      --save_reference [bool]         Save the generated reference files to the results directory
      --gencode [bool]                Use fc_group_features_type = 'gene_type' and pass '--gencode' flag to Salmon

    Strandedness:
      --forward_stranded [bool]       The library is forward stranded
      --reverse_stranded [bool]       The library is reverse stranded
      --unstranded [bool]             The default behaviour

    UMI handling:
      --with_umi [bool]               Enable UMI-tools processing steps
      --umitools_extract_method [str] The "extract method" used in the UMI tools extract step
      --umitools_bc_pattern [str]     Pattern for barcodes on read1
      --umitools_extract_extra [str]  Extra argument string which is literally passed to `umitools extract`
      --umitools_dedup_extra [str]    Extra argument string which is literally passed to `umitools dedup`
      --save_umi_intermeds [bool]     Save FastQ files with UMIs added to the read name and deduplicated BAM filesl to the results directory

    Trimming:
      --clip_r1 [int]                 Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]                 Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]     Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]     Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
      --trim_nextseq [int]            Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails
      --pico [bool]                   Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forward_stranded --clip_r1 3 --three_prime_clip_r2 3
      --skip_trimming [bool]          Skip Trim Galore step
      --save_trimmed [bool]           Save trimmed FastQ file intermediates

    Ribosomal RNA removal:
      --remove_ribo_rna [bool]        Removes ribosomal RNA using SortMeRNA
      --save_non_ribo_reads [bool]    Save FastQ file intermediates after removing rRNA
      --ribo_database_manifest [file] Path to file that contains file paths for rRNA databases, optional

    Alignment:
      --aligner [str]                 Specifies the aligner to use (available are: 'hisat2', 'star')
      --pseudo_aligner [str]          Specifies the pseudo aligner to use (available are: 'salmon'). Runs in addition to `--aligner`
      --star_align_options [str]      Additional options that will be appended to the STAR alignment command
      --hisat2_align_options [str]    Additional options that will be appended to the HISAT2 alignment command
      --stringtie_ignore_gtf [bool]   Perform reference-guided de novo assembly of transcripts using StringTie i.e. dont restrict to those in GTF file
      --seq_center [str]              Add sequencing center in @RG line of output BAM header
      --save_align_intermeds  [bool]  Save the BAM files from the aligment step - not done by default
      --save_unaligned [bool]         Save unaligned reads from either STAR, HISAT2 or Salmon to extra output files.
      --skip_rsem [bool]              Skip the RSEM step for read quantification
      --skip_alignment [bool]         Skip alignment altogether (usually in favor of pseudoalignment)
      --percent_aln_skip [float]      Percentage alignment below which samples are removed from further processing. Default: 5%

    Read counting:
      --fc_extra_attributes [str]     Define which extra parameters should also be included in featureCounts (default: 'gene_name')
      --fc_group_features [str]       Define the attribute type used to group features. (default: 'gene_id')
      --fc_count_type [str]           Define the type used to assign reads. (default: 'exon')
      --fc_group_features_type [str]  Define the type attribute used to group features based on the group attribute (default: 'gene_biotype')

    QC:
      --skip_qc [bool]                Skip all QC steps apart from MultiQC
      --skip_fastqc [bool]            Skip FastQC
      --skip_preseq [bool]            Skip Preseq
      --skip_dupradar [bool]          Skip dupRadar (and Picard MarkDuplicates)
      --skip_qualimap [bool]          Skip Qualimap
      --skip_biotype_qc [bool]        Skip Biotype QC
      --skip_rseqc [bool]             Skip RSeQC
      --skip_edger [bool]             Skip edgeR MDS plot and heatmap
      --skip_multiqc [bool]           Skip MultiQC

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false

ch_mdsplot_header = Channel.fromPath("$baseDir/assets/multiqc/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header = Channel.fromPath("$baseDir/assets/multiqc/heatmap_header.txt", checkIfExists: true)
ch_biotypes_header = Channel.fromPath("$baseDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)
Channel
    .fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
    .into { ch_where_trim_galore
            ch_where_star
            ch_where_hisat2
            ch_where_hisat2_sort
            ch_where_umi_extract
            ch_where_umi_dedup }

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
if (params.pico) {
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}

// Get rRNA databases
// Default is set to bundled DB list in `assets/rrna-db-defaults.txt`
rRNA_database = file(params.ribo_database_manifest)
if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
Channel
    .from(rRNA_database.readLines())
    .map { row -> file(row) }
    .set { sortmerna_fasta }

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2') {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if (params.pseudo_aligner && params.pseudo_aligner != 'salmon') {
    exit 1, "Invalid pseudo aligner option: ${params.pseudo_aligner}. Valid options: 'salmon'"
}

if (params.star_index && params.aligner == 'star' && !params.skip_alignment) {
    if (hasExtension(params.star_index, 'gz')) {
        Channel
            .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
            .set { star_index_gz }
    } else {
        Channel
            .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
            .set { star_index }
    }
} else if (params.hisat2_index && params.aligner == 'hisat2' && !params.skip_alignment) {
    if (hasExtension(params.hisat2_index, 'gz')) {
        Channel
            .fromPath("${params.hisat2_index}", checkIfExists: true)
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
            .set { hs2_indices_gz }
    } else {
        Channel
            .fromPath("${params.hisat2_index}*", checkIfExists: true)
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
            .set { hs2_indices }
    }
}
if (params.fasta && !params.skip_alignment) {
    if (params.additional_fasta) {
        if (hasExtension(params.additional_fasta, "gz")) {
            Channel
                .fromPath(params.additional_fasta)
                .ifEmpty { exit 1, "Additional Fasta file not found: ${params.additional_fasta}" }
                .set { additional_fasta_gz }
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
                .set { genome_fasta_gz }
        } else {
            Channel
                .fromPath(params.additional_fasta)
                .ifEmpty { exit 1, "Additional Fasta file not found: ${params.additional_fasta}" }
                .into { ch_additional_fasta_for_gtf
                        ch_additional_fasta_to_concat }
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome Fasta file not found: ${params.fasta}" }
                .set { ch_genome_fasta }
        }
    } else {
        if (hasExtension(params.fasta, "gz")) {
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
                .set { genome_fasta_gz }
        } else {
            Channel
                .fromPath(params.fasta)
                .ifEmpty { exit 1, "Genome Fasta file not found: ${params.fasta}" }
                .into { ch_fasta_for_star_index
                        ch_fasta_for_hisat_index
                        ch_fasta_for_rsem_reference }
        }
    }
} else if (params.skip_alignment) {
    println "Skipping alignment ..."
} else {
    exit 1, "No reference genome files specified!"
}

if (params.aligner == 'hisat2' && params.splicesites) {
    Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .into { indexing_splicesites
                alignment_splicesites }
}

// Separately check for whether salmon needs a genome fasta to extract
// transcripts from, or can use a transcript fasta directly
if (params.pseudo_aligner == 'salmon') {
    if (params.salmon_index) {
        if (hasExtension(params.salmon_index, 'gz')) {
            Channel
                .fromPath(params.salmon_index, checkIfExists: true)
                .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
                .set { salmon_index_gz }
        } else {
            Channel
                .fromPath(params.salmon_index, checkIfExists: true)
                .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
                .set { salmon_index }
        }
    } else if (params.transcript_fasta) {
        if (hasExtension(params.transcript_fasta, 'gz')) {
            Channel
                .fromPath(params.transcript_fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
                .set { transcript_fasta_gz }
        } else {
            Channel
                .fromPath(params.transcript_fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
                .set { ch_fasta_for_salmon_index }
        }
    } else if (params.fasta) {
        if (params.additional_fasta) {
            Channel
                .fromPath(params.additional_fasta)
                .ifEmpty { exit 1, "Additional Fasta file not found: ${params.additional_fasta}" }
                .into { ch_additional_fasta_for_gtf
                        ch_additional_fasta_to_concat }
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome Fasta file not found: ${params.fasta}" }
                .set { ch_genome_fasta }
    } else if (params.fasta && (params.gff || params.gtf)) {
        // Need to extract transcripts out of genome fasta + gtf to get
        // transcript fasta
        log.info "Extracting transcript fastas from genome fasta + gtf/gff"
        if (params.compressed_reference) {
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
                .set { genome_fasta_gz }
        } else {
            Channel
                .fromPath(params.fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
                .set { ch_fasta_for_salmon_transcripts }
        }
    }
    } else {
        exit 1, "To use with `--pseudo_aligner 'salmon'`, must provide either --transcript_fasta or both --fasta and --gtf"
    }
}

skip_rsem = params.skip_rsem
if (!params.skip_alignment && !params.skip_rsem && params.aligner != "star") {
    skip_rsem = true
    println "RSEM only works with STAR. Disabling RSEM."
}
if (params.rsem_reference && !params.skip_rsem && !params.skip_alignment) {
    Channel
        .fromPath(params.rsem_reference, checkIfExists: true)
        .ifEmpty {exit 1, "RSEM reference not found: ${params.rsem_reference}"}
        .set { rsem_reference }
}
if (params.fasta && !params.skip_alignment) {
    if (hasExtension(params.fasta, 'gz')) {
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
            .set { genome_fasta_gz }
    }
} else if (params.skip_alignment) {
    println "Skipping Alignment ..."
} else {
    exit 1, "No reference genome files specified! "
}

if (params.gtf) {
    if (params.gff) {
        // Prefer gtf over gff
        log.info "Both GTF and GFF have been provided: Using GTF as priority."
    }
    if (hasExtension(params.gtf, 'gz')) {
        Channel
            .fromPath(params.gtf, checkIfExists: true)
            .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
            .set { gtf_gz }
    } else {
        Channel
            .fromPath(params.gtf, checkIfExists: true)
            .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
            .set { gtfFile }
    }
} else if (params.gff) {
    if (hasExtension(params.gff, 'gz')) {
        Channel
            .fromPath(params.gff, checkIfExists: true)
            .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
            .set { gff_gz }
    } else {
        Channel
            .fromPath(params.gff, checkIfExists: true)
            .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
            .set { gffFile }
    }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}

if (params.bed12) {
    Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .set { bed_rseqc }
}

if (params.gencode) {
    biotype = "gene_type"
} else {
    biotype = params.fc_group_features_type
}

if (params.skip_alignment && !params.pseudo_aligner) {
    exit 1, "--skip_alignment specified without --pseudo_aligner .. did you mean to specify --pseudo_aligner salmon"
}

if (workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel') {
    if (!params.project) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (params.input_paths) {
    if (params.single_end) {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
            .into { ch_read_files_fastqc
                    raw_reads_umitools }
    } else {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
            .into { ch_read_files_fastqc
                    raw_reads_umitools }
    }
} else {
    Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_read_files_fastqc
                raw_reads_umitools }
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision)  summary['Pipeline Release'] = workflow.revision
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Input']        = params.input
summary['Data Type']    = params.single_end ? 'Single-End' : 'Paired-End'
if (params.genome)      summary['Genome'] = params.genome
if (params.pico)        summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness'] = (unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None')
summary['Trimming']     = "5'R1: $clip_r1 / 5'R2: $clip_r2 / 3'R1: $three_prime_clip_r1 / 3'R2: $three_prime_clip_r2 / NextSeq Trim: $params.trim_nextseq"
if (params.with_umi) {
    summary["With UMI"]                           = params.with_umi
    summary["umi_tools extract-method"]           = params.umitools_extract_method
    summary["umi_tools bc-pattern"]               = params.umitools_bc_pattern
    summary["umi_tools extract extra parameters"] = params.umitools_extract_extra
    summary["umi_tools dedup extra parameters"]   = params.umitools_dedup_extra
}
if (params.additional_fasta) summary["Additional Fasta"] = params.additional_fasta
if (params.aligner == 'star') {
    summary['Aligner'] = "STAR"
    if (params.star_align_options) summary['STAR Align Options'] = params.star_align_options
    if (params.star_index_options) summary['STAR Index Options'] = params.star_index_options
    if (params.star_index)         summary['STAR Index'] = params.star_index
    else if (params.fasta)         summary['Fasta Ref']  = params.fasta
} else if (params.aligner == 'hisat2') {
    summary['Aligner'] = "HISAT2"
    if (params.hisat2_align_options) summary['HISAT2 Align Options'] = params.hisat2_align_options
    if (params.hisat2_index) summary['HISAT2 Index'] = params.hisat2_index
    else if (params.fasta)   summary['Fasta Ref']    = params.fasta
    if (params.splicesites)  summary['Splice Sites'] = params.splicesites
}
if (params.pseudo_aligner == 'salmon') {
    summary['Pseudo Aligner'] = "Salmon"
    if (params.transcript_fasta) summary['Transcript Fasta'] = params.transcript_fasta
}
if (params.gtf)                    summary['GTF Annotation'] = params.gtf
if (params.gff)                    summary['GFF3 Annotation'] = params.gff
if (params.bed12)                  summary['BED Annotation'] = params.bed12
if (params.gencode)                summary['GENCODE'] = params.gencode
if (params.stringtie_ignore_gtf)   summary['StringTie Ignore GTF'] = params.stringtie_ignore_gtf
summary['Remove Ribosomal RNA']    = params.remove_ribo_rna
if (params.fc_group_features_type) summary['Biotype GTF field'] = biotype
summary['Save prefs'] = "Ref Genome: "+(params.save_reference ? 'Yes' : 'No')+" / Trimmed FastQ: "+(params.save_trimmed ? 'Yes' : 'No')+" / Alignment intermediates: "+(params.save_align_intermeds ? 'Yes' : 'No')
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']  = params.outdir
summary['Launch dir']  = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir']  = workflow.projectDir
summary['User']        = workflow.userName
if (workflow.profile == 'awsbatch') {
    summary['AWS Region']     = params.awsregion
    summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

compressed_reference = (hasExtension(params.fasta, 'gz') ||
    hasExtension(params.transcript_fasta, 'gz') || hasExtension(params.star_index, 'gz') ||
    hasExtension(params.hisat2_index, 'gz') || hasExtension(params.additional_fasta, "gz"))

if (compressed_reference) {
    // This complex logic is to prevent accessing the genome_fasta_gz variable if
    // necessary indices for STAR, HiSAT2, Salmon already exist, or if
    // params.transcript_fasta is provided as then the transcript sequences don't
    // need to be extracted.
    need_star_index = params.aligner == 'star' && !params.star_index
    need_hisat2_index = params.aligner == 'hisat2' && !params.hisat2_index
    need_rsem_ref = !params.skip_rsem && !params.rsem_reference
    //when an additional fasta is provided, the fasta and gtf file need
    //to be unzipeed to be merged in a later stage. --> Execute the following code block.
    need_aligner_index = need_hisat2_index || need_star_index || need_rsem_ref || params.additional_fasta
    alignment_no_indices = !params.skip_alignment && need_aligner_index
    pseudoalignment_no_indices = params.pseudo_aligner == "salmon" && !(params.transcript_fasta || params.salmon_index)
    if (params.fasta && (alignment_no_indices || pseudoalignment_no_indices)) {
        process GUNZIP_GENOME_FASTA {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from genome_fasta_gz

            output:
            path "${gz.baseName}" into ch_genome_fasta

            script:
            """
            gunzip --verbose --stdout --force $gz > ${gz.baseName}
            """
        }

        if (params.additional_fasta) {
            process GUNZIP_ADDITIONAL_FASTA {
                tag "$gz"
                publishDir path: { params.save_reference ? "${params.outdir}/reference/transcriptome" : params.outdir },
                    saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

                input:
                path gz from additional_fasta_gz

                output:
                path "${gz.baseName}" into ch_additional_fasta_for_gtf, ch_additional_fasta_to_concat

                script:
                """
                gunzip --verbose --stdout --force $gz > ${gz.baseName}
                """
            }
        }
    }
    if (params.gtf) {
        process GUNZIP_GTF {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from gtf_gz

            output:
            path "${gz.baseName}" into gtfFile

            script:
            """
            gunzip --verbose --stdout --force $gz > ${gz.baseName}
            """
        }
    }
    if (params.gff && !params.gtf) {
        process GUNZIP_GFF {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from gff_gz

            output:
            path "${gz.baseName}" into gffFile

            script:
            """
            gunzip --verbose --stdout --force $gz > ${gz.baseName}
            """
        }
    }
    if (params.transcript_fasta && params.pseudo_aligner == 'salmon' && !params.salmon_index) {
        process GUNZIP_TRANSCRIPT_FASTA {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/transcriptome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from transcript_fasta_gz

            output:
            path "${gz.baseName}" into ch_fasta_for_salmon_index

            script:
            """
            gunzip --verbose --stdout --force $gz > ${gz.baseName}
            """
        }
    }
    if (params.bed12) {
        process GUNZIP_BED12 {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from bed12_gz

            output:
            path "${gz.baseName}" into bed_rseqc

            script:
            """
            gunzip --verbose --stdout --force $gz > ${gz.baseName}
            """
        }
    }
    if (!params.skip_alignment && params.star_index && params.aligner == "star") {
        process GUNZIP_STAR_INDEX {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome/star" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from star_index_gz

            output:
            path "${gz.simpleName}" into star_index

            script:
            // Use tar as the star indices are a folder, not a file
            """
            tar -xzvf $gz
            """
        }
    }
    if (!params.skip_alignment && params.hisat2_index && params.aligner == 'hisat2') {
        process GUNZIP_HISAT2_INDEX {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome/hisat2" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from hs2_indices_gz

            output:
            path "*.ht2*" into hs2_indices

            script:
            // Use tar as the hisat2 indices are a folder, not a file
            """
            tar -xzvf $gz
            """
        }
    }
    if (params.salmon_index && params.pseudo_aligner == 'salmon') {
        process GUNZIP_SALMON_INDEX {
            tag "$gz"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/transcriptome/hisat2" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gz from salmon_index_gz

            output:
            path "${gz.simpleName}" into salmon_index

            script:
            // Use tar as the hisat2 indices are a folder, not a file
            """
            tar -xzvf $gz
            """
        }
    }
}

if (params.additional_fasta) {
    process MAKE_ADDITIONAL_GTF {
        tag "$fasta"
        publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path fasta from ch_additional_fasta_for_gtf

        output:
        path "${fasta.baseName}.gtf" into ch_additional_gtf

        """
        fasta2gtf.py -o ${fasta.baseName}.gtf $fasta
        """
    }

    process COMBINE_GENOME_ANNOTATIONS {
        tag "$genome_name"
        publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path genome_fasta from ch_genome_fasta
        path genome_gtf from gtfFile
        path additional_fasta from ch_additional_fasta_to_concat.collect()
        path additional_gtf from ch_additional_gtf.collect()

        output:
        path "${genome_name}.fa" into ch_fasta_for_star_index,
                                      ch_fasta_for_hisat_index,
                                      ch_fasta_for_salmon_transcripts,
                                      ch_fasta_for_rsem_reference
        path "${genome_name}.gtf" into gtf_makeSTARindex,
                                       gtf_makeHisatSplicesites,
                                       gtf_makeHISATindex,
                                       gtf_makeSalmonIndex,
                                       gtf_makeBED12,
                                       gtf_star,
                                       gtf_dupradar,
                                       gtf_featureCounts,
                                       gtf_stringtieFPKM,
                                       gtf_salmon,
                                       gtf_salmon_merge,
                                       gtf_qualimap,
                                       gtf_makeRSEMReference

        script:
        main_genome_name = params.genome ? params.genome : genome_fasta.getBaseName()
        transgenomes = additional_fasta.collect{ it.getBaseName() }.sort().join("+")
        genome_name = "${main_genome_name}_${transgenomes}"
        """
        cat $genome_fasta $additional_fasta > ${genome_name}.fa
        cat $genome_gtf $additional_gtf > ${genome_name}.gtf
        """
    }
} else {
    gtfFile
        .into { gtf_makeSTARindex
                gtf_makeHisatSplicesites
                gtf_makeHISATindex
                gtf_makeSalmonIndex
                gtf_makeBED12
                gtf_star
                gtf_dupradar
                gtf_featureCounts
                gtf_stringtieFPKM
                gtf_salmon
                gtf_salmon_merge
                gtf_qualimap
                gtf_makeRSEMReference }
}

/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if (params.gff && !params.gtf) {
    process GFF2GTF {
        tag "$gff"
        publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path gff from gffFile

        output:
        path "${gff.baseName}.gtf" into gtf_makeSTARindex,
                                        gtf_makeHisatSplicesites,
                                        gtf_makeHISATindex,
                                        gtf_makeSalmonIndex,
                                        gtf_makeBED12,
                                        gtf_star,
                                        gtf_dupradar,
                                        gtf_featureCounts,
                                        gtf_stringtieFPKM,
                                        gtf_salmon,
                                        gtf_salmon_merge,
                                        gtf_qualimap

        script:
        """
        gffread $gff --keep-exon-attrs -F -T -o ${gff.baseName}.gtf
        """
    }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if (!params.bed12) {
    process GTF2BED {
        tag "$gtf"
        publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path gtf from gtf_makeBED12

        output:
        path "${gtf.baseName}.bed" into bed_rseqc

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}

/*
 * PREPROCESSING - Build STAR index
 */
if (!params.skip_alignment) {
    if (params.aligner == 'star' && !params.star_index && params.fasta) {
        process STAR_GENOMEGENERATE {
            tag "$fasta"
            label 'high_memory'
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path fasta from ch_fasta_for_star_index
            path gtf from gtf_makeSTARindex

            output:
            path "star" into star_index

            script:
            def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
            """
            mkdir star
            STAR \\
                --runMode genomeGenerate \\
                --runThreadN $task.cpus \\
                --sjdbGTFfile $gtf \\
                --genomeDir star/ \\
                --genomeFastaFiles $fasta \\
                $avail_mem \\
                $params.star_index_options
            """
        }
    }

    /*
    * PREPROCESSING - Build HISAT2 splice sites file
    */
    if (params.aligner == 'hisat2' && !params.splicesites) {
        process MAKE_HISAT2_SPLICESITES {
            tag "$gtf"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path gtf from gtf_makeHisatSplicesites

            output:
            path "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites,
                                                                alignment_splicesites

            script:
            """
            hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
            """
        }
    }

    /*
    * PREPROCESSING - Build HISAT2 index
    */
    if (params.aligner == 'hisat2' && !params.hisat2_index && params.fasta) {
        process HISAT2_BUILD {
            tag "$fasta"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path fasta from ch_fasta_for_hisat_index
            path indexing_splicesites from indexing_splicesites
            path gtf from gtf_makeHISATindex

            output:
            path "${fasta.baseName}.*.ht2*" into hs2_indices

            script:
            if (!task.memory) {
                log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
                avail_mem = 0
            } else {
                log.info "[HISAT2 index build] Available memory: ${task.memory}"
                avail_mem = task.memory.toGiga()
            }
            if (avail_mem > params.hisat_build_memory) {
                log.info "[HISAT2 index build] Over ${params.hisat_build_memory} GB available, so using splice sites and exons in HISAT2 index"
                extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
                ss = "--ss $indexing_splicesites"
                exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
            } else {
                log.info "[HISAT2 index build] Less than ${params.hisat_build_memory} GB available, so NOT using splice sites and exons in HISAT2 index."
                log.info "[HISAT2 index build] Use --hisat_build_memory [small number] to skip this check."
                extract_exons = ''
                ss = ''
                exon = ''
            }
            """
            $extract_exons
            hisat2-build -p $task.cpus $ss $exon $fasta ${fasta.baseName}.hisat2_index
            """
        }
    }

    /**
    * PREPROCESSING - Build RSEM reference
    */
    if (!params.skip_rsem && !params.rsem_reference && params.fasta) {
        process RSEM_PREPAREREFERENCE {
            tag "$fasta"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

            input:
            path fasta from ch_fasta_for_rsem_reference
            path gtf from gtf_makeRSEMReference

            output:
            path "rsem" into rsem_reference

            script:
            """
            mkdir rsem
            rsem-prepare-reference -p $task.cpus --gtf $gtf $fasta rsem/ref
            """
        }
    }
}


/*
 * PREPROCESSING - Create Salmon transcriptome index
 */
if (params.pseudo_aligner == 'salmon' && !params.salmon_index) {
    if (!params.transcript_fasta) {
        process TRANSCRIPTS_TO_FASTA {
            tag "$fasta"
            publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
                saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode


            input:
            path fasta from ch_fasta_for_salmon_transcripts
            path gtf from gtf_makeSalmonIndex

            output:
            path "*.fa" into ch_fasta_for_salmon_index

            script:
	          // filter_gtf_for_genes_in_genome.py is bundled in this package, in rnaseq/bin
            """
            filter_gtf_for_genes_in_genome.py --gtf $gtf --fasta $fasta -o ${gtf.baseName}__in__${fasta.baseName}.gtf
            gffread -F -w transcripts.fa -g $fasta ${gtf.baseName}__in__${fasta.baseName}.gtf
            """
        }
    }
    process SALMON_INDEX {
        tag "$fasta"
        label "salmon"
        publishDir path: { params.save_reference ? "${params.outdir}/reference/genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path fasta from ch_fasta_for_salmon_index

        output:
        path 'salmon_index' into salmon_index

        script:
        def gencode = params.gencode  ? "--gencode" : ""
        """
        salmon index --threads $task.cpus -t $fasta $gencode -i salmon_index
        """
    }
}

/*
 * STEP 1 - FastQC
 */
process FASTQC {
    tag "$name"
    label 'mid_memory'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
            filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
        }

    when:
    !params.skip_fastqc

    input:
    tuple val(name), path(reads) from ch_read_files_fastqc

    output:
    path "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    def threads = params.single_end ? 1 : 2
    """
    fastqc --quiet --threads $threads $reads
    """
}

/*
 * STEP 1+ - UMItools
 */
if (params.with_umi) {
    process UMITOOLS_EXTRACT {
        tag "$name"
        label "low_memory"
        publishDir "${params.outdir}/umitools/extract", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.endsWith('.log')) filename
                else if (!params.save_umi_intermeds && filename == "where_are_my_files.txt") filename
                else if (params.save_umi_intermeds && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        tuple val(name), path(reads) from raw_reads_umitools
        path wherearemyfiles from ch_where_umi_extract.collect()

        output:
        tuple val(name), path("*fq.gz") into raw_reads_trimgalore
        path "*.log"
        path "where_are_my_files.txt"

        script:
        if (params.single_end) {
            """
            umi_tools extract \\
                -I $reads \\
                -S ${name}_umi_extracted.fq.gz \\
                --extract-method=${params.umitools_extract_method} \\
                --bc-pattern="${params.umitools_bc_pattern}" \\
                ${params.umitools_extract_extra} > ${name}_umi_extract.log
            """
        }  else {
            """
            umi_tools extract \\
                -I ${reads[0]} \\
                --read2-in=${reads[1]} \\
                -S ${name}_umi_extracted_R1.fq.gz \\
                --read2-out=${name}_umi_extracted_R2.fq.gz \\
                --extract-method=${params.umitools_extract_method} \\
                --bc-pattern="${params.umitools_bc_pattern}" \\
                ${params.umitools_extract_extra} > ${name}_umi_extract.log
            """
        }
    }
} else {
    raw_reads_trimgalore = raw_reads_umitools
    umi_tools_extract_results = Channel.empty()
}


/*
 * STEP 2 - Trim Galore!
 */
if (!params.skip_trimming) {
    process TRIMGALORE {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/trimgalore", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("_fastqc") > 0) "fastqc/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else if (!params.save_trimmed && filename == "where_are_my_files.txt") filename
                else if (params.save_trimmed && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        tuple val(name), path(reads) from raw_reads_trimgalore
        path wherearemyfiles from ch_where_trim_galore.collect()

        output:
        tuple val(name), path("*fq.gz") into trimgalore_reads
        path "*trimming_report.txt" into trimgalore_results
        path "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
        path "where_are_my_files.txt"

        script:
        // Calculate number of --cores for TrimGalore based on value of task.cpus
        // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
        // See: https://github.com/nf-core/atacseq/pull/65
        def cores = 1
        if (task.cpus) {
            cores = (task.cpus as int) - 4
            if (params.single_end) cores = (task.cpus as int) - 3
            if (cores < 1) cores = 1
            if (cores > 4) cores = 4
        }

        c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
        c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
        tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
        tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
        nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
        if (params.single_end) {
            """
            trim_galore --cores $cores --fastqc --gzip $c_r1 $tpc_r1 $nextseq $reads
            """
        } else {
            """
            trim_galore --cores $cores --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq $reads
            """
        }
    }
} else {
    trimgalore_reads = raw_reads_trimgalore
    trimgalore_results = Channel.empty()
}

/*
 * STEP 2+ - SortMeRNA - remove rRNA sequences on request
 */
if (!params.remove_ribo_rna) {
    trimgalore_reads
        .into { trimmed_reads_alignment
                trimmed_reads_salmon }
    sortmerna_logs = Channel.empty()
} else {
    process SORTMERNA_INDEXDBRNA {
        tag "${fasta.baseName}"
        label 'low_memory'

        input:
        path fasta from sortmerna_fasta

        output:
        val "${fasta.baseName}" into sortmerna_db_name
        path "$fasta" into sortmerna_db_fasta
        path "${fasta.baseName}*" into sortmerna_db

        script:
        """
        indexdb_rna --ref $fasta,${fasta.baseName} -m 3072 -v
        """
    }

    process SORTMERNA {
        tag "$name"
        label 'low_memory'
        publishDir "${params.outdir}/sortmerna", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("_rRNA_report.txt") > 0) "logs/$filename"
                else if (params.saveNonRiboRNAReads) "reads/$filename"
                else null
            }

        input:
        tuple val(name), path(reads) from trimgalore_reads
        val db_name from sortmerna_db_name.collect()
        path db_fasta from sortmerna_db_fasta.collect()
        path db from sortmerna_db.collect()

        output:
        tuple val(name), path("*.fq.gz") into trimmed_reads_alignment,
                                              trimmed_reads_salmon
        path "*_rRNA_report.txt" into sortmerna_logs


        script:
        //concatenate reference files: ${db_fasta},${db_name}:${db_fasta},${db_name}:...
        def Refs = ''
        for (i=0; i<db_fasta.size(); i++) { Refs+= ":${db_fasta[i]},${db_name[i]}" }
        Refs = Refs.substring(1)
        if (params.single_end) {
            """
            gzip -d --force < $reads > all-reads.fastq

            sortmerna --ref $Refs \\
                --reads all-reads.fastq \\
                --num_alignments 1 \\
                -a $task.cpus \\
                --fastx \\
                --aligned rRNA-reads \\
                --other non-rRNA-reads \\
                --log -v

            gzip --force < non-rRNA-reads.fastq > ${name}.fq.gz

            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        } else {
            """
            gzip -d --force < ${reads[0]} > reads-fw.fq
            gzip -d --force < ${reads[1]} > reads-rv.fq
            merge-paired-reads.sh reads-fw.fq reads-rv.fq all-reads.fastq

            sortmerna --ref $Refs \\
                --reads all-reads.fastq \\
                --num_alignments 1 \\
                -a $task.cpus \\
                --fastx --paired_in \\
                --aligned rRNA-reads \\
                --other non-rRNA-reads \\
                --log -v

            unmerge-paired-reads.sh non-rRNA-reads.fastq non-rRNA-reads-fw.fq non-rRNA-reads-rv.fq
            gzip < non-rRNA-reads-fw.fq > ${name}-fw.fq.gz
            gzip < non-rRNA-reads-rv.fq > ${name}-rv.fq.gz

            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        }
    }
}

/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
good_alignment_scores = [:]
poor_alignment_scores = [:]
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if (percent_aligned.toFloat() <= params.percent_aln_skip.toFloat()) {
        log.info "#${c_red}################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<${c_reset}"
        poor_alignment_scores[logname] = percent_aligned
        return false
    } else {
        log.info "-${c_green}           Passed alignment > star ($logname)   >> ${percent_aligned}% <<${c_reset}"
        good_alignment_scores[logname] = percent_aligned
        return true
    }
}
if (!params.skip_alignment) {
    if (params.aligner == 'star') {
        hisat_stdout = Channel.empty()
        process STAR_ALIGN {
            tag "$name"
            label 'high_memory'
            publishDir "${params.outdir}/star", mode: params.publish_dir_mode,
                saveAs: { filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else if (params.save_unaligned && filename != "where_are_my_files.txt" && 'Unmapped' in filename) unmapped/filename
                    else if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
                    else if (params.save_align_intermeds && filename != "where_are_my_files.txt") filename
                    else null
                }

            input:
            tuple val(name), path(reads) from trimmed_reads_alignment
            path index from star_index.collect()
            path gtf from gtf_star.collect()
            path wherearemyfiles from ch_where_star.collect()

            output:
            tuple path("*Log.final.out"), path('*.sortedByCoord.out.bam'), path('*.toTranscriptome.out.bam') into star_aligned
            path "*.out" into alignment_logs
            path "*SJ.out.tab"
            path "*Log.out" into star_log
            path "where_are_my_files.txt"
            path "*Unmapped*" optional true
            path "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index

            script:
            prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
            seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
            def star_mem = task.memory ?: params.star_memory ?: false
            def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
            def unaligned = params.save_unaligned ? "--outReadsUnmapped Fastx" : ''
            """
            STAR --genomeDir $index \\
                --sjdbGTFfile $gtf \\
                --readFilesIn $reads  \\
                --runThreadN $task.cpus \\
                --twopassMode Basic \\
                --outWigType bedGraph \\
                --outSAMtype BAM SortedByCoordinate $avail_mem \\
                --readFilesCommand zcat \\
                --runDirPerm All_RWX $unaligned \\
                --quantMode TranscriptomeSAM \\
                --outFileNamePrefix $prefix $seq_center \\
                --runRNGseed 0 \\
                $params.star_align_options

            samtools index ${prefix}Aligned.sortedByCoord.out.bam
            """
        }
        // Filter removes all 'aligned' channels that fail the check
        star_bams = Channel.create()
        star_bams_transcriptome = Channel.create()
        star_aligned
            .filter { logs, bams, bams_transcriptome -> check_log(logs) }
            .separate (star_bams, star_bams_transcriptome) {
                bam_set -> [bam_set[1], bam_set[2]]
            }
        bam = star_bams
        bam_transcriptome = star_bams_transcriptome
    }

    /*
    * STEP 3 - align with HISAT2
    */
    if (params.aligner == 'hisat2') {
        star_log = Channel.empty()
        process HISAT2_ALIGN {
            tag "$name"
            label 'high_memory'
            publishDir "${params.outdir}/hisat2", mode: params.publish_dir_mode,
                saveAs: { filename ->
                    if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                    else if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
                    else if (params.save_align_intermeds && filename != "where_are_my_files.txt") filename
                    else null
                }

            input:
            tuple val(name), path(reads) from trimmed_reads_alignment
            path hs2_indices from hs2_indices.collect()
            path alignment_splicesites from alignment_splicesites.collect()
            path wherearemyfiles from ch_where_hisat2.collect()

            output:
            path "${prefix}.bam" into hisat2_bam
            path "${prefix}.hisat2_summary.txt" into alignment_logs
            path "where_are_my_files.txt"
            path "unmapped.hisat2*" optional true

            script:
            index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/
            prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
            seq_center = params.seq_center ? "--rg-id ${prefix} --rg CN:${params.seq_center.replaceAll('\\s','_')} SM:$prefix" : "--rg-id ${prefix} --rg SM:$prefix"
            def rnastrandness = ''
            if (forward_stranded && !unstranded) {
                rnastrandness = params.single_end ? '--rna-strandness F' : '--rna-strandness FR'
            } else if (reverse_stranded && !unstranded) {
                rnastrandness = params.single_end ? '--rna-strandness R' : '--rna-strandness RF'
            }
            if (params.single_end) {
                unaligned = params.save_unaligned ? "--un-gz unmapped.hisat2.gz" : ''
                """
                hisat2 \\
                    -x $index_base \\
                    -U $reads \\
                    $rnastrandness \\
                    --known-splicesite-infile $alignment_splicesites \\
                    -p $task.cpus $unaligned \\
                    --met-stderr \\
                    --new-summary \\
                    --dta \\
                    $params.hisat2_align_options \\
                    --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
                    | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
                """
            } else {
                unaligned = params.save_unaligned ? "--un-conc-gz unmapped.hisat2.gz" : ''
                """
                hisat2 \\
                    -x $index_base \\
                    -1 ${reads[0]} \\
                    -2 ${reads[1]} \\
                    $rnastrandness \\
                    --known-splicesite-infile $alignment_splicesites \\
                    --no-mixed \\
                    --no-discordant \\
                    -p $task.cpus $unaligned \\
                    --met-stderr \\
                    --new-summary \\
                    --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
                    | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
                """
            }
        }

        process HISAT2_SORT_BAM {
            tag "${bam.baseName}"
            label 'mid_memory'
            publishDir "${params.outdir}/hisat2", mode: params.publish_dir_mode,
                saveAs: { filename ->
                    if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
                    else if (params.save_align_intermeds && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
                    else null
                }

            input:
            path bam from hisat2_bam
            path wherearemyfiles from ch_where_hisat2_sort.collect()

            output:
            path "${bam.baseName}.sorted.bam" into bam
            path "${bam.baseName}.sorted.bam.bai" into bam_index
            path "where_are_my_files.txt"

            script:
            def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
            def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
            """
            samtools sort \\
                $bam \\
                -@ $task.cpus $avail_mem \\
                -o ${bam.baseName}.sorted.bam
            samtools index ${bam.baseName}.sorted.bam
            """
        }
    }

    /*
    * Step 3+ - Deduplicate bam files based on UMIs
    */
    if (params.with_umi) {
        // preseq does not work on deduplicated BAM file. Pass it the raw BAM file.
        bam
            .into { bam_umitools_dedup
                    bam_preseq }
        bam_index_umitools_dedup = bam_index

        process UMITOOLS_DEDUP {
            tag "${bam.baseName}"
            label "mid_memory"
            publishDir "${params.outdir}/umitools/dedup", mode: params.publish_dir_mode,
                saveAs: { filename ->
                    if (filename.endsWith('.tsv')) filename
                    else if (!params.save_umi_intermeds && filename == "where_are_my_files.txt") filename
                    else if (params.save_umi_intermeds && filename != "where_are_my_files.txt") filename
                    else null
                }

            input:
            path bam from bam_umitools_dedup
            path bai from bam_index_umitools_dedup
            path wherearemyfiles from ch_where_umi_dedup.collect()

            output:
            path "*.bam" into bam_dedup
            path "*.bai" into bam_dedup_index
            path "where_are_my_files.txt"
            path "*.tsv"

            script:
            """
            umi_tools dedup \\
                -I $bam \\
                -S ${bam.baseName}_deduplicated.bam \\
                --output-stats=${bam.baseName} \\
                $params.umitools_dedup_extra
            samtools index ${bam.baseName}_deduplicated.bam
            """
        }

        // RSEM transcriptome BAM file treated separately...
        if (!skip_rsem) {
            process UMITOOLS_DEDUP_TRANSCRIPTOME {
                tag "${bam.baseName}"
                label "mid_memory"
                publishDir "${params.outdir}/umitools/dedup/transcriptome", mode: params.publish_dir_mode,
                    saveAs: { filename ->
                        if (filename.endsWith('.tsv')) filename
                        else if (params.save_umi_intermeds) filename
                        else null
                    }

                input:
                path bam from bam_transcriptome

                output:
                path "*_deduplicated.bam" into bam_rsem
                path "*.tsv"

                script:
                // the transcriptome BAM file is not sorted or indexed by STAR
                // since this is the only process consuming this BAM file,
                // sorting and indexing happens right here.
                def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
                def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
                """
                samtools sort \\
                    $bam \\
                    -@ $task.cpus $avail_mem \\
                    -o ${bam.baseName}.sorted.bam
                samtools index ${bam.baseName}.sorted.bam

                umi_tools dedup \\
                    -I ${bam.baseName}.sorted.bam \\
                    -S ${bam.baseName}_deduplicated.bam \\
                    --output-stats=${bam.baseName} \\
                    $params.umitools_dedup_extra
                """
            }
        }

        bam_dedup
            .into { bam_count
                    bam_rseqc
                    bam_qualimap
                    bam_markduplicates
                    bam_featurecounts
                    bam_stringtieFPKM
                    bam_forSubsamp
                    bam_skipSubsamp }
        bam_dedup_index
            .into { bam_index_rseqc
                    bam_index_genebody }
    } else {
        bam
            .into { bam_count
                    bam_rseqc
                    bam_qualimap
                    bam_preseq
                    bam_markduplicates
                    bam_featurecounts
                    bam_stringtieFPKM
                    bam_forSubsamp
                    bam_skipSubsamp }
        bam_index
            .into { bam_index_rseqc
                    bam_index_genebody }
        if (!skip_rsem) {
            bam_rsem = bam_transcriptome
        }
    }

    /*
    * STEP 4 - RSeQC analysis
    */
    process RSEQC {
        tag "${bam.baseName - '.sorted'}"
        label 'mid_memory'
        publishDir "${params.outdir}/rseqc" , mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("bam_stat.txt") > 0)                           "bam_stat/$filename"
                else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
                else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
                else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
                else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
                else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
                else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
                else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
                else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
                else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
                else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
                else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
                else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
                else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
                else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
                else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
                else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
                else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
                else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
                else if (filename.indexOf("junction_annotation_log.txt") > 0)       "junction_annotation/$filename"
                else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
                else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
                else filename
            }

        when:
        !params.skip_qc && !params.skip_rseqc

        input:
        path bam from bam_rseqc
        path bai from bam_index_rseqc
        path bed12 from bed_rseqc.collect()

        output:
        path "*.{txt,pdf,r,xls}" into rseqc_results

        script:
        """
        infer_experiment.py -i $bam -r $bed12 > ${bam.baseName}.infer_experiment.txt
        junction_annotation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12 2> ${bam.baseName}.junction_annotation_log.txt
        bam_stat.py -i $bam 2> ${bam.baseName}.bam_stat.txt
        junction_saturation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
        inner_distance.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
        read_distribution.py -i $bam -r $bed12 > ${bam.baseName}.read_distribution.txt
        read_duplication.py -i $bam -o ${bam.baseName}.read_duplication
        """
    }

    /*
    * STEP 5 - preseq analysis
    */
    process PRESEQ {
        tag "${bam.baseName - '.sorted'}"
        label 'high_time'
        publishDir "${params.outdir}/preseq", mode: params.publish_dir_mode

        when:
        !params.skip_qc && !params.skip_preseq

        input:
        path bam from bam_preseq

        output:
        path "${bam.baseName}.ccurve.txt" into preseq_results

        script:
        """
        preseq lc_extrap -v -B $bam -o ${bam.baseName}.ccurve.txt
        """
    }

    /*
    * STEP 6 - Mark duplicates
    */
    process PICARD_MARKDUPLICATES {
        tag "${bam.baseName - '.sorted'}"
        publishDir "${params.outdir}/markduplicates", mode: params.publish_dir_mode,
            saveAs: { filename ->
                filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"
            }

        when:
        !params.skip_qc && !params.skip_dupradar

        input:
        path bam from bam_markduplicates

        output:
        path "${bam.baseName}.markDups.bam" into bam_md
        path "${bam.baseName}.markDups_metrics.txt" into picard_results
        path "${bam.baseName}.markDups.bam.bai"

        script:
        markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2)+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
        """
        picard $markdup_java_options MarkDuplicates \\
            INPUT=$bam \\
            OUTPUT=${bam.baseName}.markDups.bam \\
            METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
            REMOVE_DUPLICATES=false \\
            ASSUME_SORTED=true \\
            PROGRAM_RECORD_ID='null' \\
            VALIDATION_STRINGENCY=LENIENT
        samtools index ${bam.baseName}.markDups.bam
        """
    }

    /*
    * STEP 7 - Qualimap
    */
    process QUALIMAP {
        tag "${bam.baseName}"
        label 'low_memory'
        publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode

        when:
        !params.skip_qc && !params.skip_qualimap

        input:
        path bam from bam_qualimap
        path gtf from gtf_qualimap.collect()

        output:
        path "${bam.baseName}" into qualimap_results

        script:
        def qualimap_direction = 'non-strand-specific'
        if (forward_stranded) {
            qualimap_direction = 'strand-specific-forward'
        } else if (reverse_stranded) {
            qualimap_direction = 'strand-specific-reverse'
        }
        def paired = params.single_end ? '' : '-pe'
        def memory = task.memory.toGiga() + "G"
        """
        unset DISPLAY
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        qualimap --java-mem-size=$memory rnaseq -p $qualimap_direction $paired -bam $bam -gtf $gtf -outdir ${bam.baseName}
        """
    }

    /*
    * STEP 8 - dupRadar
    */
    process DUPRADAR {
        tag "${bam.baseName - '.sorted.markDups'}"
        label 'high_time'
        publishDir "${params.outdir}/dupradar", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
                else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
                else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
                else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
                else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
                else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
                else "$filename"
            }

        when:
        !params.skip_qc && !params.skip_dupradar

        input:
        path bam from bam_md
        path gtf from gtf_dupradar.collect()

        output:
        path "*.{pdf,txt}" into dupradar_results

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        def dupradar_direction = 0
        if (forward_stranded && !unstranded) {
            dupradar_direction = 1
        } else if (reverse_stranded && !unstranded) {
            dupradar_direction = 2
        }
        def paired = params.single_end ? 'single' :  'paired'
        """
        dupRadar.r $bam $gtf $dupradar_direction $paired $task.cpus
        """
    }

    /*
    * STEP 9 - Feature counts
    */
    process SUBREAD_FEATURECOUNTS {
        tag "${bam.baseName - '.sorted'}"
        label 'low_memory'
        publishDir "${params.outdir}/featurecounts", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
                else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
                else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
                else "$filename"
            }

        input:
        path bam from bam_featurecounts
        path gtf from gtf_featureCounts.collect()
        path biotypes_header from ch_biotypes_header.collect()

        output:
        path "${bam.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
        path "${bam.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
        path "${bam.baseName}_biotype_counts*mqc.{txt,tsv}" optional true into featureCounts_biotype

        script:
        def featureCounts_direction = 0
        def extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
        if (forward_stranded && !unstranded) {
            featureCounts_direction = 1
        } else if (reverse_stranded && !unstranded) {
            featureCounts_direction = 2
        }
        // Try to get real sample name
        sample_name = bam.baseName - 'Aligned.sortedByCoord.out' - '_subsamp.sorted'
        biotype_qc = params.skip_biotype_qc ? '' : "featureCounts -a $gtf -g $biotype -t ${params.fc_count_type} -o ${bam.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam"
        mod_biotype = params.skip_biotype_qc ? '' : "cut -f 1,7 ${bam.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam.baseName}_biotype_counts_mqc.txt && mqc_features_stat.py ${bam.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam.baseName}_biotype_counts_gs_mqc.tsv"
        """
        featureCounts \\
            -a $gtf \\
            -g $params.fc_group_features \\
            -t $params.fc_count_type \\
            -o ${bam.baseName}_gene.featureCounts.txt \\
            $extraAttributes \\
            -p \\
            -s $featureCounts_direction \\
            $bam
        $biotype_qc
        $mod_biotype
        """
    }

    /*
    * STEP 10 - Merge featurecounts
    */
    process MERGE_FEATURECOUNTS {
        tag "${input_files[0].baseName - '.sorted'}"
        label "mid_memory"
        publishDir "${params.outdir}/featurecounts", mode: params.publish_dir_mode

        input:
        path input_files from featureCounts_to_merge.collect()

        output:
        path 'merged_gene_counts.txt' into featurecounts_merged

        script:
        // Redirection (the `<()`) for the win!
        // Geneid in 1st column and gene_name in 7th
        gene_ids = "<(tail -n +2 ${input_files[0]} | cut -f1,7 )"
        counts = input_files.collect{filename ->
            // Remove first line and take third column
            "<(tail -n +2 ${filename} | sed 's:.sorted.bam::' | cut -f8)"}.join(" ")
        """
        paste $gene_ids $counts > merged_gene_counts.txt
        """
    }

    if (!skip_rsem) {
        /**
        * Step 11 - RSEM
        */
        process RSEM_CALCULATEEXPRESSION {
            tag "${bam.baseName - '.sorted'}"
            label "mid_memory"
            publishDir "${params.outdir}/rsem", mode: params.publish_dir_mode

            input:
            path bam from bam_rsem
            path "rsem" from rsem_reference.collect()

            output:
            path "*.genes.results" into rsem_results_genes
            path "*.isoforms.results" into rsem_results_isoforms
            path "*.stat" into rsem_logs

            script:
            def sample_name = bam.baseName - 'Aligned.toTranscriptome.out' - '_subsamp'
            def paired_end_flag = params.single_end ? "" : "--paired-end"
            """
            REF_FILENAME=\$(basename rsem/*.grp)
            REF_NAME="\${REF_FILENAME%.*}"
            rsem-calculate-expression \\
                -p $task.cpus \\
                $paired_end_flag \\
                --bam \\
                --estimate-rspd \\
                --append-names \\
                $bam \\
                rsem/\$REF_NAME \\
                $sample_name
            """
        }

        /**
        * Step 12 - merge RSEM TPM and counts
        */
        process MERGE_RSEM_COUNTS {
            tag "${rsem_res_gene[0].baseName}"
            label "low_memory"
            publishDir "${params.outdir}/rsem", mode: params.publish_dir_mode

            input:
            path rsem_res_gene from rsem_results_genes.collect()
            path rsem_res_isoform from rsem_results_isoforms.collect()

            output:
            path "rsem_tpm_gene.txt"
            path "rsem_tpm_isoform.txt"
            path "rsem_transcript_counts_gene.txt"
            path "rsem_transcript_counts_isoform.txt"

            script:
            """
            echo "gene_id\tgene_symbol" > gene_ids.txt
            echo "transcript_id\tgene_symbol" > transcript_ids.txt
            cut -f 1 ${rsem_res_gene[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> gene_ids.txt
            cut -f 1 ${rsem_res_isoform[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> transcript_ids.txt
            mkdir tmp_genes tmp_isoforms
            for fileid in $rsem_res_gene; do
                basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.tpm.txt
                grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_genes/\${fileid}.tpm.txt
                basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.counts.txt
                grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_genes/\${fileid}.counts.txt
            done
            for fileid in $rsem_res_isoform; do
                basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.tpm.txt
                grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_isoforms/\${fileid}.tpm.txt
                basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.counts.txt
                grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_isoforms/\${fileid}.counts.txt
            done
            paste gene_ids.txt tmp_genes/*.tpm.txt > rsem_tpm_gene.txt
            paste gene_ids.txt tmp_genes/*.counts.txt > rsem_transcript_counts_gene.txt
            paste transcript_ids.txt tmp_isoforms/*.tpm.txt > rsem_tpm_isoform.txt
            paste transcript_ids.txt tmp_isoforms/*.counts.txt > rsem_transcript_counts_isoform.txt
            """
        }
    } else {
        rsem_logs = Channel.empty()
    }

    /*
    * STEP 13 - stringtie FPKM
    */
    process STRINGTIE {
        tag "${bam.baseName - '.sorted'}"
        publishDir "${params.outdir}/stringtie", mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
                else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
                else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
                else "$filename"
            }

        input:
        path bam from bam_stringtieFPKM
        path gtf from gtf_stringtieFPKM.collect()

        output:
        path "${bam.baseName}_transcripts.gtf"
        path "${bam.baseName}.gene_abund.txt"
        path "${bam}.cov_refs.gtf"
        path "${bam.baseName}_ballgown"

        script:
        def st_direction = ''
        if (forward_stranded && !unstranded) {
            st_direction = "--fr"
        } else if (reverse_stranded && !unstranded) {
            st_direction = "--rf"
        }
        def ignore_gtf = params.stringtie_ignore_gtf ? "" : "-e"
        """
        stringtie $bam \\
            $st_direction \\
            -o ${bam.baseName}_transcripts.gtf \\
            -v \\
            -G $gtf \\
            -A ${bam.baseName}.gene_abund.txt \\
            -C ${bam}.cov_refs.gtf \\
            -b ${bam.baseName}_ballgown \\
            $ignore_gtf
        """
    }

    /*
    * STEP 14 - edgeR MDS and heatmap
    */
    process SAMPLE_CORRELATION {
        tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
        label 'low_memory'
        publishDir "${params.outdir}/sample_correlation", mode: params.publish_dir_mode

        when:
        !params.skip_qc && !params.skip_edger

        input:
        path input_files from geneCounts.collect()
        val num_bams from bam_count.count()
        path mdsplot_header from ch_mdsplot_header
        path heatmap_header from ch_heatmap_header

        output:
        path "*.{txt,pdf,csv}" into sample_correlation_results

        when:
        num_bams > 2 && (!params.sample_level)

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        edgeR_heatmap_MDS.r $input_files
        cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
        mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
        cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
        mv tmp_file log2CPM_sample_correlation_mqc.csv
        """
    }
} else {
    star_log = Channel.empty()
    hisat_stdout = Channel.empty()
    alignment_logs = Channel.empty()
    rseqc_results = Channel.empty()
    picard_results = Channel.empty()
    qualimap_results = Channel.empty()
    sample_correlation_results = Channel.empty()
    featureCounts_logs = Channel.empty()
    dupradar_results = Channel.empty()
    preseq_results = Channel.empty()
    featureCounts_biotype = Channel.empty()
    rsem_logs = Channel.empty()
}


/*
 * STEP 15 - Transcriptome quantification with Salmon
 */
if (params.pseudo_aligner == 'salmon') {
    process SALMON_QUANT {
        tag "$sample"
        publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode

        input:
        tuple val(sample), path(reads) from trimmed_reads_salmon
        path index from salmon_index.collect()
        path gtf from gtf_salmon.collect()

        output:
        path "${sample}/" into salmon_logs
        tuple val(sample), path("${sample}/") into salmon_tximport,
                                                   salmon_parsegtf

        script:
        def rnastrandness = params.single_end ? 'U' : 'IU'
        if (forward_stranded && !unstranded) {
            rnastrandness = params.single_end ? 'SF' : 'ISF'
        } else if (reverse_stranded && !unstranded) {
            rnastrandness = params.single_end ? 'SR' : 'ISR'
        }
        def endedness = params.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
        def unmapped = params.save_unaligned ? "--writeUnmappedNames" : ''
        """
        salmon quant \\
            --validateMappings \\
            --seqBias --useVBOpt --gcBias \\
            --geneMap ${gtf} \\
            --threads $task.cpus \\
            --libType=${rnastrandness} \\
            --index ${index} \\
            $endedness $unmapped\\
            -o ${sample}
        """
    }

    process SALMON_TX2GENE {
        label 'low_memory'
        publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode

        input:
        path ("salmon/*") from salmon_parsegtf.collect()
        path gtf from gtf_salmon_merge

        output:
        path "tx2gene.csv" into salmon_tx2gene, salmon_merge_tx2gene

        script:
        """
        parse_gtf.py --gtf $gtf --salmon salmon --id $params.fc_group_features --extra $params.fc_extra_attributes -o tx2gene.csv
        """
    }

    process SALMON_TXIMPORT {
        label 'low_memory'
        publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode

        input:
        tuple val(name), path("salmon/*") from salmon_tximport
        path tx2gene from salmon_tx2gene.collect()

        output:
        path "${name}_salmon_gene_tpm.csv" into salmon_gene_tpm
        path "${name}_salmon_gene_counts.csv" into salmon_gene_counts
        path "${name}_salmon_transcript_tpm.csv" into salmon_transcript_tpm
        path "${name}_salmon_transcript_counts.csv" into salmon_transcript_counts

        script:
        """
        tximport.r NULL salmon $name
        """
    }

    process SALMON_MERGE {
        label 'mid_memory'
        publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode

        input:
        path gene_tpm_files from salmon_gene_tpm.collect()
        path gene_count_files from salmon_gene_counts.collect()
        path transcript_tpm_files from salmon_transcript_tpm.collect()
        path transcript_count_files from salmon_transcript_counts.collect()
        path tx2gene from salmon_merge_tx2gene

        output:
        path "salmon_merged*.csv" into salmon_merged_ch
        path "*.rds"

        script:
        // First field is the gene/transcript ID
        gene_ids = "<(cut -f1 -d, ${gene_tpm_files[0]} | tail -n +2 | cat <(echo '${params.fc_group_features}') - )"
        transcript_ids = "<(cut -f1 -d, ${transcript_tpm_files[0]} | tail -n +2 | cat <(echo 'transcript_id') - )"

        // Second field is counts/TPM
        gene_tpm = gene_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
        gene_counts = gene_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
        transcript_tpm = transcript_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
        transcript_counts = transcript_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
        """
        paste -d, $gene_ids $gene_tpm > salmon_merged_gene_tpm.csv
        paste -d, $gene_ids $gene_counts > salmon_merged_gene_counts.csv
        paste -d, $transcript_ids $transcript_tpm > salmon_merged_transcript_tpm.csv
        paste -d, $transcript_ids $transcript_counts > salmon_merged_transcript_counts.csv

        se.r NULL salmon_merged_gene_counts.csv salmon_merged_gene_tpm.csv
        se.r NULL salmon_merged_transcript_counts.csv salmon_merged_transcript_tpm.csv
        """
    }
} else {
    salmon_logs = Channel.empty()
}

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/rnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    path 'software_versions_mqc.yaml' into ch_software_versions_yaml
    path "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    sortmerna --version &> v_sortmerna.txt
    STAR --version &> v_star.txt
    hisat2 --version &> v_hisat2.txt
    stringtie --version &> v_stringtie.txt
    preseq &> v_preseq.txt
    read_duplication.py --version &> v_rseqc.txt
    bamCoverage --version &> v_deeptools.txt || true
    featureCounts -v &> v_featurecounts.txt
    rsem-calculate-expression --version &> v_rsem.txt
    salmon --version &> v_salmon.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
    Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    umi_tools --version &> v_umi_tools.txt
    unset DISPLAY && qualimap rnaseq &> v_qualimap.txt || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 16 - MultiQC
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    when:
    !params.skip_multiqc

    input:
    path multiqc_config from ch_multiqc_config
    path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    path ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    path ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    path ('alignment/*') from alignment_logs.collect().ifEmpty([])
    path ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    path ('picard/*') from picard_results.collect().ifEmpty([])
    path ('qualimap/*') from qualimap_results.collect().ifEmpty([])
    path ('preseq/*') from preseq_results.collect().ifEmpty([])
    path ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    path ('featurecounts/*') from featureCounts_logs.collect().ifEmpty([])
    path ('featurecounts_biotype/*') from featureCounts_biotype.collect().ifEmpty([])
    path ('rsem/*') from rsem_logs.collect().ifEmpty([])
    path ('salmon/*') from salmon_logs.collect().ifEmpty([])
    path ('sample_correlation/*') from sample_correlation_results.collect().ifEmpty([])
    path ('sortmerna/*') from sortmerna_logs.collect().ifEmpty([])
    path ('software_versions/*') from ch_software_versions_yaml.collect()
    path workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    path "*multiqc_report.html" into ch_multiqc_report
    path "*_data"
    path "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}

/*
 * STEP 17 - Output Description HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path output_docs from ch_output_docs
    path images from ch_output_docs_images

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[nf-core/rnaseq] Successful: $workflow.runName"
    if (poor_alignment_scores.size() > 0) {
        subject = "[nf-core/rnaseq] Partially Successful (${poor_alignment_scores.size()} skipped): $workflow.runName"
    }
    if (!workflow.success) {
      subject = "[nf-core/rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = poor_alignment_scores.keySet()
    email_fields['percent_aln_skip'] = params.percent_aln_skip
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skip_multiqc) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/rnaseq] Found multiple reports from process 'MULTIQC', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/rnaseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/rnaseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if (mqc_report.size() <= params.max_multiqc_email_size.toBytes()) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/rnaseq] Sent summary e-mail to $email_address (mail)"
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

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (good_alignment_scores.size() > 0) {
        total_aln_count = good_alignment_scores.size() + poor_alignment_scores.size()
        idx = 0;
        samp_aln = ''
        for (samp in good_alignment_scores) {
            samp_aln += "    ${samp.key}: ${samp.value}%\n"
            idx += 1
            if (idx > 5) {
                samp_aln += "    ..see pipeline reports for full list\n"
                break;
            }
        }
        log.info "[${c_purple}nf-core/rnaseq${c_reset}] ${c_green}${good_alignment_scores.size()}/$total_aln_count samples passed minimum ${params.percent_aln_skip}% aligned check\n${samp_aln}${c_reset}"
    }
    if (poor_alignment_scores.size() > 0) {
        samp_aln = ''
        poor_alignment_scores.each { samp, value ->
            samp_aln += "    ${samp}: ${value}%\n"
        }
        log.info "[${c_purple}nf-core/rnaseq${c_reset}] ${c_red} WARNING - ${poor_alignment_scores.size()} samples skipped due to poor mapping percentages!\n${samp_aln}${c_reset}"
    }

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "- ${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "- ${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "- ${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "[${c_purple}nf-core/rnaseq${c_reset}] ${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "[${c_purple}nf-core/rnaseq${c_reset}] ${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnaseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
