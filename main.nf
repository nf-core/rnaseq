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

nextflow.enable.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

/*
 * Reference genomes
 */
// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.fasta        = params.genome ? params.genomes[ params.genome ].fasta  ?: false : false
params.gtf          = params.genome ? params.genomes[ params.genome ].gtf    ?: false : false
params.gff          = params.genome ? params.genomes[ params.genome ].gff    ?: false : false
params.bed12        = params.genome ? params.genomes[ params.genome ].bed12  ?: false : false
params.star_index   = params.genome ? params.genomes[ params.genome ].star   ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
anno_readme         = params.genome ? params.genomes[ params.genome ].readme ?: false : false

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check mandatory parameters
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, 'Genome fasta file not specified!' }
if (!params.gtf && !params.gff) { exit 1, "No GTF or GFF3 annotation specified!" }
if (params.gtf && params.gff)   { log.info "Both GTF and GFF have been provided: Using GTF as priority." }

// Check input path parameters to see if they exist
checkPathParamList = [
    params.gtf, params.gff, params.bed12, params.additional_fasta, params.transcript_fasta,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index,
    params.splicesites, params.ribo_database_manifest
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check rRNA databases for sortmerna
ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
if (ch_ribo_db.isEmpty()) {exit 1, "File ${ch_ribo_db.getName()} is empty!"}

// Check alignment parameters
if (params.skip_alignment) { log.info "Skipping alignment processes..." }
if (params.aligner != 'star' && params.aligner != 'hisat2') {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if (params.pseudo_aligner && params.pseudo_aligner != 'salmon') {
    exit 1, "Invalid pseudo aligner option: ${params.pseudo_aligner}. Valid options: 'salmon'"
}
if (params.skip_alignment && !params.pseudo_aligner) {
    exit 1, "--skip_alignment specified without --pseudo_aligner .. did you mean to specify --pseudo_aligner salmon"
}
// if (params.pseudo_aligner == 'salmon') {
//     if (!params.salmon_index || !params.transcript_fasta || (!params.fasta || !(params.gff || params.gtf))) {
//         exit 1, "To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff"
//     }
// }

// Check if we are running RSEM
def skip_rsem = params.skip_rsem
if (!params.skip_rsem) {
    if (params.aligner != "star") {
        skip_rsem = true
        log.info "RSEM only works when '--aligner star' is set. Disabling RSEM."
    }
}

// Save AWS IGenomes file containing annotation version
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
 * Initiate parameters
 */
// Define regular variables so that they can be overwritten
def forward_stranded = params.forward_stranded
def reverse_stranded = params.reverse_stranded
def unstranded = params.unstranded
if (params.pico) {
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}

// Set biotype for featureCounts
def biotype = params.gencode ? "gene_type" : params.fc_group_features_type

/*
 * Check other parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

// if (workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel') {
//     if (!params.project) exit 1, "No UPPMAX project ID found! Use --project"
// }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images    = file("$baseDir/docs/images/", checkIfExists: true)

// Header files for MultiQC
ch_mdsplot_header        = file("$baseDir/assets/multiqc/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header        = file("$baseDir/assets/multiqc/heatmap_header.txt", checkIfExists: true)
ch_biotypes_header       = file("$baseDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}

summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { UNTAR as UNTAR_RSEM_INDEX                      } from './modules/local/process/untar'
include { RSEM_PREPAREREFERENCE                          } from './modules/local/process/rsem_preparereference'
include { CAT_FASTQ                                      } from './modules/local/process/cat_fastq'
include { SORTMERNA                                      } from './modules/local/process/sortmerna'
include { UMITOOLS_DEDUP as UMITOOLS_DEDUP_GENOME
          UMITOOLS_DEDUP as UMITOOLS_DEDUP_TRANSCRIPTOME } from './modules/local/process/umitools_dedup'
include { OUTPUT_DOCUMENTATION                           } from './modules/local/process/output_documentation'
include { GET_SOFTWARE_VERSIONS                          } from './modules/local/process/get_software_versions'
// include { MULTIQC                                        } from './modules/local/process/multiqc'

include { INPUT_CHECK                                    } from './modules/local/subworkflow/input_check'
include { PREP_GENOME                                    } from './modules/local/subworkflow/prep_genome'
include { ALIGN_STAR                                     } from './modules/local/subworkflow/align_star'
include { ALIGN_HISAT2                                   } from './modules/local/subworkflow/align_hisat2'
include { QUANTIFY_SALMON                                } from './modules/local/subworkflow/quantify_salmon'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { SAMTOOLS_INDEX             } from './modules/nf-core/software/samtools/index/main'
include { PRESEQ_LCEXTRAP            } from './modules/nf-core/software/preseq/lcextrap/main'
// include { SUBREAD_FEATURECOUNTS      } from './modules/nf-core/software/subread/featurecounts/main'

include { FASTQC_UMITOOLS_TRIMGALORE } from './modules/nf-core/subworkflow/fastqc_umitools_trimgalore'
include { BAM_SORT_SAMTOOLS          } from './modules/nf-core/subworkflow/bam_sort_samtools'
include { MARK_DUPLICATES_PICARD     } from './modules/nf-core/subworkflow/mark_duplicates_picard'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
def good_alignment_scores = [:]
def poor_alignment_scores = [:]
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - '.Log.final'
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

workflow {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    def publish_genome_options = params.save_reference ? [publish_dir : 'genome']       : [publish_files : [:]]
    def publish_index_options  = params.save_reference ? [publish_dir : 'genome/index'] : [publish_files : [:]]
    PREP_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.bed12,
        params.additional_fasta,
        publish_genome_options
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREP_GENOME.out.gffread_version.ifEmpty(null))

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input, params.seq_center, [:] )
        .map {
            meta, bam ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_cat_fastq }

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ ( ch_cat_fastq, params.modules['cat_fastq'] )

    /*
     * SUBWORKFLOW: Read QC, extract UMI and trim adapters
     */
    def method = params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : ''
    def pattern = params.umitools_bc_pattern ? "--bc-pattern='${params.umitools_bc_pattern}'" : ''
    params.modules['umitools_extract'].args += " $method $pattern"
    if (params.save_umi_intermeds) { params.modules['umitools_extract'].publish_files.put('fastq.gz','') }

    def nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    if (params.save_trimmed) { params.modules['trimgalore'].publish_files.put('fq.gz','') }

    FASTQC_UMITOOLS_TRIMGALORE (
        CAT_FASTQ.out.reads,
        params.skip_fastqc,
        params.with_umi,
        params.skip_trimming,
        params.modules['fastqc'],
        params.modules['umitools_extract'],
        params.modules['trimgalore']
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.umitools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    /*
     * MODULE: Remove ribosomal RNA reads
     */
    ch_trimmed_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads
    ch_sortmerna_log = Channel.empty()
    if (params.remove_ribo_rna) {
        if (params.save_non_ribo_reads) { params.modules['sortmerna'].publish_files.put('fastq.gz','') }
        ch_sortmerna_fasta = Channel.from(ch_ribo_db.readLines()).map { row -> file(row) }.collect()
        SORTMERNA ( ch_trimmed_reads, ch_sortmerna_fasta, params.modules['sortmerna'] )
            .reads
            .set { ch_trimmed_reads }
        ch_sortmerna_log = SORTMERNA.out.log
        ch_software_versions = ch_software_versions.mix(SORTMERNA.out.version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Alignment with STAR
     */
    ch_genome_bam        = Channel.empty()
    ch_genome_bai        = Channel.empty()
    ch_transcriptome_bam = Channel.empty()
    ch_star_log          = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star') {
        // TODO nf-core: Not working - only save indices if --save_reference is specified
        if (params.save_reference)       { params.modules['star_genomegenerate']['publish_files'] = null }
        if (params.save_align_intermeds) {
            params.modules['star_align'].publish_files.put('bam','')
            params.modules['samtools_sort'].publish_dir += '/star'
            params.modules['samtools_sort'].publish_files = ['bam':'', 'bai':'', 'stats':'samtools_stats', 'flagstat':'samtools_stats', 'idxstats':'samtools_stats']
        }
        if (params.save_unaligned)       { params.modules['star_align'].publish_files.put('fastq.gz','unmapped') }
        def unaligned = params.save_unaligned ? " --outReadsUnmapped Fastx" : ''
        params.modules['star_align'].args += unaligned

        ALIGN_STAR (
            ch_trimmed_reads,
            params.star_index,
            PREP_GENOME.out.fasta,
            PREP_GENOME.out.gtf,
            params.modules['star_genomegenerate'],
            params.modules['star_align'],
            params.modules['samtools_sort']
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bai        = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_star_log          = ALIGN_STAR.out.log_final
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))

        //     // Filter removes all 'aligned' channels that fail the check
        //     star_bams = Channel.create()
        //     star_bams_transcriptome = Channel.create()
        //     star_aligned
        //         .filter { logs, bams, bams_transcriptome -> check_log(logs) }
        //         .separate (star_bams, star_bams_transcriptome) {
        //             bam_set -> [bam_set[1], bam_set[2]]
        //         }
        //     bam = star_bams
        //     bam_transcriptome = star_bams_transcriptome
        // }
        //
        // /*
        //  * SUBWORKFLOW: Gene/transcript quantification with RSEM
        //  */
        // if (!skip_rsem) {
        //     if (params.rsem_index) {
        //         if (params.rsem_index.endsWith('.tar.gz')) {
        //             ch_rsem_index = UNTAR_RSEM_INDEX ( params.rsem_index, publish_index_options ).untar
        //         } else {
        //             ch_rsem_index = file(params.rsem_index)
        //         }
        //     } else {
        //         // TODO nf-core: Not working - only save indices if --save_reference is specified
        //         if (params.save_reference) { params.modules['rsem_preparereference']['publish_files'] = null }
        //         //println(params.modules['rsem_preparereference'])
        //         ch_rsem_index = RSEM_PREPAREREFERENCE ( ch_fasta, ch_gtf, params.modules['rsem_preparereference'] ).index
        //     }
        // }
    }

    /*
     * SUBWORKFLOW: Alignment with HISAT2
     */
    ch_hisat2_log = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2' && params.aligner != 'star') {
        // TODO nf-core: Not working - only save indices if --save_reference is specified
        if (params.save_reference)       { params.modules['hisat2_build']['publish_files'] = null }
        if (params.save_align_intermeds) {
            params.modules['hisat2_align'].publish_files.put('bam','')
            params.modules['samtools_sort'].publish_dir += '/hisat2'
            params.modules['samtools_sort'].publish_files = ['bam':'', 'bai':'', 'stats':'samtools_stats', 'flagstat':'samtools_stats', 'idxstats':'samtools_stats']
        }
        if (params.save_unaligned)       { params.modules['hisat2_align'].publish_files.put('fastq.gz','unmapped') }

        ALIGN_HISAT2 (
            ch_trimmed_reads,
            params.hisat2_index,
            PREP_GENOME.out.fasta,
            PREP_GENOME.out.gtf,
            params.splicesites,
            params.modules['hisat2_build'],
            params.modules['hisat2_align'],
            params.modules['samtools_sort']
        )
        ch_genome_bam = ALIGN_HISAT2.out.bam
        ch_genome_bai = ALIGN_HISAT2.out.bai
        ch_hisat2_log = ALIGN_HISAT2.out.summary
        ch_software_versions = ch_software_versions.mix(ALIGN_HISAT2.out.hisat2_version.first().ifEmpty(null))
    }

    /*
     * MODULE: Run Preseq
     */
    if (!params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP ( ch_genome_bam, params.modules['preseq_lcextrap'] )
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))
        //preseq lc_extrap -v -B $bam -o ${bam.baseName}.ccurve.txt
    }
    //
    //     process PICARD_MARKDUPLICATES {
    //         tag "${bam.baseName - '.sorted'}"
    //         publishDir "${params.outdir}/markduplicates", mode: params.publish_dir_mode,
    //             saveAs: { filename ->
    //                 filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"
    //             }
    //
    //         when:
    //         !params.skip_qc && !params.skip_dupradar
    //
    //         input:
    //         path bam from bam_markduplicates
    //
    //         output:
    //         path "${bam.baseName}.markDups.bam" into bam_md
    //         path "${bam.baseName}.markDups_metrics.txt" into picard_results
    //         path "${bam.baseName}.markDups.bam.bai"
    //
    //         script:
    //         markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2)+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
    //         """
    //         picard $markdup_java_options MarkDuplicates \\
    //             INPUT=$bam \\
    //             OUTPUT=${bam.baseName}.markDups.bam \\
    //             TMP_DIR='./tmp' \\
    //             METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
    //             REMOVE_DUPLICATES=false \\
    //             ASSUME_SORTED=true \\
    //             PROGRAM_RECORD_ID='null' \\
    //             VALIDATION_STRINGENCY=LENIENT
    //         samtools index ${bam.baseName}.markDups.bam
    //         """
    //     }

    /*
     * MODULE: Remove duplicate reads based on UMIs
     */
    if (!params.skip_alignment && params.with_umi) {
        if (params.save_umi_intermeds) {
            params.modules['umitools_dedup_genome'].publish_files.put('bam','')
            params.modules['umitools_dedup_genome'].publish_files.put('bai','')
        }
        UMITOOLS_DEDUP_GENOME ( ch_genome_bam.join(ch_genome_bai, by: [0]), params.modules['umitools_dedup_genome'] )
        SAMTOOLS_INDEX ( UMITOOLS_DEDUP_GENOME.out.bam, params.modules['umitools_dedup_genome'] )
        ch_genome_bam = UMITOOLS_DEDUP_GENOME.out.bam
        ch_genome_bai = SAMTOOLS_INDEX.out.bai

        if (!skip_rsem) {
            BAM_SORT_SAMTOOLS ( ch_transcriptome_bam, params.modules['samtools_sort_umitools_dedup'] )

            if (params.save_umi_intermeds) { params.modules['umitools_dedup_transcriptome'].publish_files.put('bam','') }
            UMITOOLS_DEDUP_TRANSCRIPTOME (
                BAM_SORT_SAMTOOLS.out.bam.join(BAM_SORT_SAMTOOLS.out.bai, by: [0]),
                params.modules['umitools_dedup_transcriptome']
            )
            ch_transcriptome_bam = UMITOOLS_DEDUP_TRANSCRIPTOME.out.bam
        }
    }

    /*
     * SUBWORKFLOW: Mark duplicate reads
     */
    // if (!params.skip_alignment && !params.skip_markduplicates) {
    //     // MARK_DUPLICATES_PICARD (
    //     //     ch_genome_bam,
    //     //     params.modules['picard_markduplicates'],
    //     //     params.modules['picard_markduplicates_samtools']
    //     // )
    //     // ch_software_versions = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    //     //
    //     //         picard $markdup_java_options MarkDuplicates \\
    //     //             INPUT=$bam \\
    //     //             OUTPUT=${bam.baseName}.markDups.bam \\
    //     //             TMP_DIR='./tmp' \\
    //     //             METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
    //     //             REMOVE_DUPLICATES=false \\
    //     //             ASSUME_SORTED=true \\
    //     //             PROGRAM_RECORD_ID='null' \\
    //     //             VALIDATION_STRINGENCY=LENIENT
    //     //         samtools index ${bam.baseName}.markDups.bam
    // }

    /*
     * SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
     */
    ch_salmon_log = Channel.empty()
    if (params.pseudo_aligner == 'salmon') {
        // TODO nf-core: Not working - only save indices if --save_reference is specified
        if (params.save_reference) { params.modules['salmon_index']['publish_files'] = null }
        def gencode = params.gencode  ? " --gencode" : ""
        params.modules['salmon_index'].args += gencode

        def unmapped = params.save_unaligned ? " --writeUnmappedNames" : ''
        params.modules['salmon_quant'].args += unmapped

        QUANTIFY_SALMON (
            ch_trimmed_reads,
            params.salmon_index,
            PREP_GENOME.out.fasta,
            PREP_GENOME.out.gtf,
            publish_index_options,
            publish_genome_options,
            params.modules['salmon_index'],
            params.modules['salmon_quant']
        )
        ch_salmon_log = QUANTIFY_SALMON.out.results
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.version.first().ifEmpty(null))
    }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect(), [publish_files : ['csv':'']] )
    OUTPUT_DOCUMENTATION  ( ch_output_docs, ch_output_docs_images, [:] )
//
//     /*
//      * MultiQC
//      */
//     workflow_summary = Schema.params_mqc_summary(summary)
//     ch_workflow_summary = Channel.value(workflow_summary)
//     MULTIQC (
//         ch_multiqc_config,
//         ch_multiqc_custom_config.collect().ifEmpty([]),
//         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
//         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
//
//         FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
//         FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
//         FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
//
//         MAP_BWA_MEM.out.stats.collect{it[1]},
//         MAP_BWA_MEM.out.flagstat.collect{it[1]},
//         MAP_BWA_MEM.out.idxstats.collect{it[1]},
//
//         MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
//         MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
//         MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
//         MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),
//
//         PRESEQ_LCEXTRAP.out.ccurve.collect{it[1]}.ifEmpty([]),
//
//         SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]),
//
//         params.modules['multiqc']
//     )

}

//
// ////////////////////////////////////////////////////
// /* --              COMPLETION EMAIL            -- */
// ////////////////////////////////////////////////////
//
// // Find a way to pass MultiQC report here as well as poor_alignment_scores, good_alignment_scores
// workflow.onComplete {
//     def multiqc_report = []
//     Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log, poor_alignment_scores)
//     Completion.summary(workflow, params, log, poor_alignment_scores, good_alignment_scores)
// }
//
////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
//
//     if (params.with_umi) {//
//         bam_dedup
//             .into { bam_count
//                     bam_rseqc
//                     bam_qualimap
//                     bam_markduplicates
//                     bam_featurecounts
//                     bam_stringtieFPKM
//                     bam_forSubsamp
//                     bam_skipSubsamp }
//         bam_dedup_index
//             .into { bam_index_rseqc
//                     bam_index_genebody }
//     } else {
//         bam
//             .into { bam_count
//                     bam_rseqc
//                     bam_qualimap
//                     bam_preseq
//                     bam_markduplicates
//                     bam_featurecounts
//                     bam_stringtieFPKM
//                     bam_forSubsamp
//                     bam_skipSubsamp }
//         bam_index
//             .into { bam_index_rseqc
//                     bam_index_genebody }
//         if (!skip_rsem) {
//             bam_rsem = bam_transcriptome
//         }
//     }
//
//     process RSEQC {
//         tag "${bam.baseName - '.sorted'}"
//         label 'mid_memory'
//         publishDir "${params.outdir}/rseqc" , mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("bam_stat.txt") > 0)                           "bam_stat/$filename"
//                 else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
//                 else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
//                 else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
//                 else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
//                 else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
//                 else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
//                 else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
//                 else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
//                 else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
//                 else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
//                 else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
//                 else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
//                 else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
//                 else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
//                 else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
//                 else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
//                 else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
//                 else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
//                 else if (filename.indexOf("junction_annotation_log.txt") > 0)       "junction_annotation/$filename"
//                 else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
//                 else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
//                 else filename
//             }
//
//         when:
//         !params.skip_qc && !params.skip_rseqc
//
//         input:
//         path bam from bam_rseqc
//         path bai from bam_index_rseqc
//         path bed12 from ch_bed12
//
//         output:
//         path "*.{txt,pdf,r,xls}" into rseqc_results
//
//         script:
//         """
//         infer_experiment.py -i $bam -r $bed12 > ${bam.baseName}.infer_experiment.txt
//         junction_annotation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12 2> ${bam.baseName}.junction_annotation_log.txt
//         bam_stat.py -i $bam 2> ${bam.baseName}.bam_stat.txt
//         junction_saturation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
//         inner_distance.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
//         read_distribution.py -i $bam -r $bed12 > ${bam.baseName}.read_distribution.txt
//         read_duplication.py -i $bam -o ${bam.baseName}.read_duplication
//         """
//     }
//
//     process QUALIMAP {
//         tag "${bam.baseName}"
//         label 'low_memory'
//         publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode
//
//         when:
//         !params.skip_qc && !params.skip_qualimap
//
//         input:
//         path bam from bam_qualimap
//         path gtf from ch_gtf
//
//         output:
//         path "${bam.baseName}" into qualimap_results
//
//         script:
//         def qualimap_direction = 'non-strand-specific'
//         if (forward_stranded) {
//             qualimap_direction = 'strand-specific-forward'
//         } else if (reverse_stranded) {
//             qualimap_direction = 'strand-specific-reverse'
//         }
//         def paired = params.single_end ? '' : '-pe'
//         def memory = task.memory.toGiga() + "G"
//         """
//         unset DISPLAY
//         export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
//         qualimap --java-mem-size=$memory rnaseq -p $qualimap_direction $paired -bam $bam -gtf $gtf -outdir ${bam.baseName}
//         """
//     }
//
//     process DUPRADAR {
//         tag "${bam.baseName - '.sorted.markDups'}"
//         label 'high_time'
//         publishDir "${params.outdir}/dupradar", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
//                 else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
//                 else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
//                 else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
//                 else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
//                 else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
//                 else "$filename"
//             }
//
//         when:
//         !params.skip_qc && !params.skip_dupradar
//
//         input:
//         path bam from bam_md
//         path gtf from ch_gtf
//
//         output:
//         path "*.{pdf,txt}" into dupradar_results
//
//         script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
//         def dupradar_direction = 0
//         if (forward_stranded && !unstranded) {
//             dupradar_direction = 1
//         } else if (reverse_stranded && !unstranded) {
//             dupradar_direction = 2
//         }
//         def paired = params.single_end ? 'single' :  'paired'
//         """
//         dupRadar.r $bam $gtf $dupradar_direction $paired $task.cpus
//         """
//     }
//
//     process SUBREAD_FEATURECOUNTS {
//         tag "${bam.baseName - '.sorted'}"
//         label 'low_memory'
//         publishDir "${params.outdir}/featurecounts", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
//                 else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
//                 else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
//                 else "$filename"
//             }
//
//         input:
//         path bam from bam_featurecounts
//         path gtf from ch_gtf
//         path biotypes_header from ch_biotypes_header
//
//         output:
//         path "${bam.baseName}_gene.featureCounts.txt" into geneCounts,
//                                                            featureCounts_to_merge
//         path "${bam.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
//         path "${bam.baseName}_biotype_counts*mqc.{txt,tsv}" optional true into featureCounts_biotype
//
//         script:
//         def featureCounts_direction = 0
//         def extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
//         if (forward_stranded && !unstranded) {
//             featureCounts_direction = 1
//         } else if (reverse_stranded && !unstranded) {
//             featureCounts_direction = 2
//         }
//         // Try to get real sample name
//         sample_name = bam.baseName - 'Aligned.sortedByCoord.out' - '_subsamp.sorted'
//         biotype_qc = params.skip_biotype_qc ? '' : "featureCounts -a $gtf -g $biotype -t ${params.fc_count_type} -o ${bam.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam"
//         mod_biotype = params.skip_biotype_qc ? '' : "cut -f 1,7 ${bam.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam.baseName}_biotype_counts_mqc.txt && mqc_features_stat.py ${bam.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam.baseName}_biotype_counts_gs_mqc.tsv"
//         """
//         featureCounts \\
//             -a $gtf \\
//             -g $params.fc_group_features \\
//             -t $params.fc_count_type \\
//             -o ${bam.baseName}_gene.featureCounts.txt \\
//             $extraAttributes \\
//             -p \\
//             -s $featureCounts_direction \\
//             $bam
//         $biotype_qc
//         $mod_biotype
//         """
//     }
//
//     process MERGE_FEATURECOUNTS {
//         tag "${input_files[0].baseName - '.sorted'}"
//         label "mid_memory"
//         publishDir "${params.outdir}/featurecounts", mode: params.publish_dir_mode
//
//         input:
//         path input_files from featureCounts_to_merge.collect()
//
//         output:
//         path 'merged_gene_counts.txt' into featurecounts_merged
//
//         script:
//         // Redirection (the `<()`) for the win!
//         // Geneid in 1st column and gene_name in 7th
//         gene_ids = "<(tail -n +2 ${input_files[0]} | cut -f1,7 )"
//         counts = input_files.collect{filename ->
//             // Remove first line and take third column
//             "<(tail -n +2 ${filename} | sed 's:.sorted.bam::' | cut -f8)"}.join(" ")
//         """
//         paste $gene_ids $counts > merged_gene_counts.txt
//         """
//     }
//
//     if (!skip_rsem) {
//         process RSEM_CALCULATEEXPRESSION {
//             tag "${bam.baseName - '.sorted'}"
//             label "mid_memory"
//             publishDir "${params.outdir}/rsem", mode: params.publish_dir_mode
//
//             input:
//             path bam from bam_rsem
//             path "rsem" from ch_rsem_index
//
//             output:
//             path "*.genes.results" into rsem_results_genes
//             path "*.isoforms.results" into rsem_results_isoforms
//             path "*.stat" into rsem_logs
//
//             script:
//             def sample_name = bam.baseName - '.Aligned.toTranscriptome.out' - '_subsamp'
//             def paired_end = params.single_end ? "" : "--paired-end"
//             """
//             REF_FILENAME=\$(basename rsem/*.grp)
//             REF_NAME="\${REF_FILENAME%.*}"
//             rsem-calculate-expression \\
//                 -p $task.cpus \\
//                 $paired_end \\
//                 --bam \\
//                 --estimate-rspd \\
//                 --append-names \\
//                 $bam \\
//                 rsem/\$REF_NAME \\
//                 $sample_name
//             """
//         }
//
//         process MERGE_RSEM_COUNTS {
//             tag "${rsem_res_gene[0].baseName}"
//             label "low_memory"
//             publishDir "${params.outdir}/rsem", mode: params.publish_dir_mode
//
//             input:
//             path rsem_res_gene from rsem_results_genes.collect()
//             path rsem_res_isoform from rsem_results_isoforms.collect()
//
//             output:
//             path "rsem_tpm_gene.txt"
//             path "rsem_tpm_isoform.txt"
//             path "rsem_transcript_counts_gene.txt"
//             path "rsem_transcript_counts_isoform.txt"
//
//             script:
//             """
//             echo "gene_id\tgene_symbol" > gene_ids.txt
//             echo "transcript_id\tgene_symbol" > transcript_ids.txt
//             cut -f 1 ${rsem_res_gene[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> gene_ids.txt
//             cut -f 1 ${rsem_res_isoform[0]} | grep -v "^#" | tail -n+2 | sed -E "s/(_PAR_Y)?(_|\$)/\\1\\t/" >> transcript_ids.txt
//             mkdir tmp_genes tmp_isoforms
//             for fileid in $rsem_res_gene; do
//                 basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.tpm.txt
//                 grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_genes/\${fileid}.tpm.txt
//                 basename \$fileid | sed s/\\.genes.results\$//g > tmp_genes/\${fileid}.counts.txt
//                 grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_genes/\${fileid}.counts.txt
//             done
//             for fileid in $rsem_res_isoform; do
//                 basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.tpm.txt
//                 grep -v "^#" \${fileid} | cut -f 6 | tail -n+2 >> tmp_isoforms/\${fileid}.tpm.txt
//                 basename \$fileid | sed s/\\.isoforms.results\$//g > tmp_isoforms/\${fileid}.counts.txt
//                 grep -v "^#" \${fileid} | cut -f 5 | tail -n+2 >> tmp_isoforms/\${fileid}.counts.txt
//             done
//             paste gene_ids.txt tmp_genes/*.tpm.txt > rsem_tpm_gene.txt
//             paste gene_ids.txt tmp_genes/*.counts.txt > rsem_transcript_counts_gene.txt
//             paste transcript_ids.txt tmp_isoforms/*.tpm.txt > rsem_tpm_isoform.txt
//             paste transcript_ids.txt tmp_isoforms/*.counts.txt > rsem_transcript_counts_isoform.txt
//             """
//         }
//     } else {
//         rsem_logs = Channel.empty()
//     }
//
//     process STRINGTIE {
//         tag "${bam.baseName - '.sorted'}"
//         publishDir "${params.outdir}/stringtie", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
//                 else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
//                 else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
//                 else "$filename"
//             }
//
//         input:
//         path bam from bam_stringtieFPKM
//         path gtf from ch_gtf
//
//         output:
//         path "${bam.baseName}_transcripts.gtf"
//         path "${bam.baseName}.gene_abund.txt"
//         path "${bam}.cov_refs.gtf"
//         path "${bam.baseName}_ballgown"
//
//         script:
//         def st_direction = ''
//         if (forward_stranded && !unstranded) {
//             st_direction = "--fr"
//         } else if (reverse_stranded && !unstranded) {
//             st_direction = "--rf"
//         }
//         def ignore_gtf = params.stringtie_ignore_gtf ? "" : "-e"
//         """
//         stringtie $bam \\
//             $st_direction \\
//             -o ${bam.baseName}_transcripts.gtf \\
//             -v \\
//             -G $gtf \\
//             -A ${bam.baseName}.gene_abund.txt \\
//             -C ${bam}.cov_refs.gtf \\
//             -b ${bam.baseName}_ballgown \\
//             $ignore_gtf
//         """
//     }
//
//     process SAMPLE_CORRELATION {
//         tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - '.Aligned'}"
//         label 'low_memory'
//         publishDir "${params.outdir}/sample_correlation", mode: params.publish_dir_mode
//
//         when:
//         !params.skip_qc && !params.skip_edger
//
//         input:
//         path input_files from geneCounts.collect()
//         val num_bams from bam_count.count()
//         path mdsplot_header from ch_mdsplot_header
//         path heatmap_header from ch_heatmap_header
//
//         output:
//         path "*.{txt,pdf,csv}" into sample_correlation_results
//
//         when:
//         num_bams > 2 && (!params.sample_level)
//
//         script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
//         """
//         edgeR_heatmap_MDS.r $input_files
//         cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
//         mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
//         cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
//         mv tmp_file log2CPM_sample_correlation_mqc.csv
//         """
//     }
// } else {
//     star_log = Channel.empty()
//     hisat_stdout = Channel.empty()
//     alignment_logs = Channel.empty()
//     rseqc_results = Channel.empty()
//     picard_results = Channel.empty()
//     qualimap_results = Channel.empty()
//     sample_correlation_results = Channel.empty()
//     featureCounts_logs = Channel.empty()
//     dupradar_results = Channel.empty()
//     preseq_results = Channel.empty()
//     featureCounts_biotype = Channel.empty()
//     rsem_logs = Channel.empty()
// }
//
// process GET_SOFTWARE_VERSIONS {
//     publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//             if (filename.indexOf(".csv") > 0) filename
//             else null
//         }
//
//     output:
//     path 'software_versions_mqc.yaml' into ch_software_versions_yaml
//     path "software_versions.csv"
//
//     script:
//     """
//     echo $workflow.manifest.version &> v_ngi_rnaseq.txt
//     echo $workflow.nextflow.version &> v_nextflow.txt
//     fastqc --version &> v_fastqc.txt
//     cutadapt --version &> v_cutadapt.txt
//     trim_galore --version &> v_trim_galore.txt
//     sortmerna --version &> v_sortmerna.txt
//     STAR --version &> v_star.txt
//     hisat2 --version &> v_hisat2.txt
//     stringtie --version &> v_stringtie.txt
//     preseq &> v_preseq.txt
//     read_duplication.py --version &> v_rseqc.txt
//     bamCoverage --version &> v_deeptools.txt || true
//     featureCounts -v &> v_featurecounts.txt
//     rsem-calculate-expression --version &> v_rsem.txt
//     salmon --version &> v_salmon.txt
//     picard MarkDuplicates --version &> v_markduplicates.txt  || true
//     samtools --version &> v_samtools.txt
//     multiqc --version &> v_multiqc.txt
//     Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
//     Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
//     umi_tools --version &> v_umi_tools.txt
//     unset DISPLAY && qualimap rnaseq &> v_qualimap.txt || true
//     scrape_software_versions.py &> software_versions_mqc.yaml
//     """
// }
//
// process MULTIQC {
//     publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode
//
//     when:
//     !params.skip_multiqc
//
//     input:
//     path multiqc_config from ch_multiqc_config
//     path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
//     path ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
//     path ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
//     path ('alignment/*') from alignment_logs.collect().ifEmpty([])
//     path ('rseqc/*') from rseqc_results.collect().ifEmpty([])
//     path ('picard/*') from picard_results.collect().ifEmpty([])
//     path ('qualimap/*') from qualimap_results.collect().ifEmpty([])
//     path ('preseq/*') from preseq_results.collect().ifEmpty([])
//     path ('dupradar/*') from dupradar_results.collect().ifEmpty([])
//     path ('featurecounts/*') from featureCounts_logs.collect().ifEmpty([])
//     path ('featurecounts_biotype/*') from featureCounts_biotype.collect().ifEmpty([])
//     path ('rsem/*') from rsem_logs.collect().ifEmpty([])
//     path ('salmon/*') from salmon_logs.collect().ifEmpty([])
//     path ('sample_correlation/*') from sample_correlation_results.collect().ifEmpty([])
//     path ('sortmerna/*') from sortmerna_logs.collect().ifEmpty([])
//     path ('software_versions/*') from ch_software_versions_yaml.collect()
//     path workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
//
//     output:
//     path "*multiqc_report.html" into ch_multiqc_report
//     path "*_data"
//     path "multiqc_plots"
//
//     script:
//     rtitle = run_name ? "--title \"$run_name\"" : ''
//     rfilename = run_name ? "--filename " + run_name.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     """
//     multiqc . -f $rtitle $rfilename $custom_config_file
//     """
// }
