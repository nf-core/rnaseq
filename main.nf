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

ch_where_are_my_files    = file("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)

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

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_BED12
    GUNZIP as GUNZIP_ADDITIONAL_FASTA
    GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from './modules/local/process/gunzip'
include {
    UNTAR as UNTAR_STAR_INDEX
    UNTAR as UNTAR_HISAT2_INDEX
    UNTAR as UNTAR_RSEM_INDEX
    UNTAR as UNTAR_SALMON_INDEX       } from './modules/local/process/untar'
include { GFFREAD                     } from './modules/local/process/gffread'
include { GTF2BED                     } from './modules/local/process/gtf2bed'
include { CAT_ADDITIONAL_FASTA        } from './modules/local/process/cat_additional_fasta'
include { STAR_GENOMEGENERATE         } from './modules/local/process/star_genomegenerate'
include { HISAT2_EXTRACTSPLICESITES   } from './modules/local/process/hisat2_extractsplicesites'
include { HISAT2_BUILD                } from './modules/local/process/hisat2_build'
include { RSEM_PREPAREREFERENCE       } from './modules/local/process/rsem_preparereference'
include { TRANSCRIPTS2FASTA           } from './modules/local/process/transcripts2fasta'
include { SALMON_INDEX                } from './modules/local/process/salmon_index'
//include { GET_CHROM_SIZES             } from './modules/local/process/get_chrom_sizes'
include { OUTPUT_DOCUMENTATION        } from './modules/local/process/output_documentation'
include { GET_SOFTWARE_VERSIONS       } from './modules/local/process/get_software_versions'
// include { MULTIQC                     } from './modules/local/process/multiqc'

include { INPUT_CHECK                 } from './modules/local/subworkflow/input_check'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

//include { PICARD_COLLECTMULTIPLEMETRICS } from './modules/nf-core/software/picard/collectmultiplemetrics/main'
//include { PRESEQ_LCEXTRAP               } from './modules/nf-core/software/preseq/lcextrap/main'
//include { UCSC_BEDRAPHTOBIGWIG          } from './modules/nf-core/software/ucsc/bedgraphtobigwig/main'
//include { SUBREAD_FEATURECOUNTS         } from './modules/nf-core/software/subread/featurecounts/main'

include { FASTQC_TRIMGALORE             } from './modules/nf-core/subworkflow/fastqc_trimgalore'
//include { MARK_DUPLICATES_PICARD        } from './modules/nf-core/subworkflow/mark_duplicates_picard'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /*
     * PREPROCESSING: Initialise variables / empty channels
     */
    ch_software_versions = Channel.empty()

    /*
     * PREPROCESSING: Uncompress genome fasta file if required
     */
    def publish_genome = params.save_reference ? [publish_dir : 'genome'] : [publish_files : [:]]
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta, publish_genome ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * PREPROCESSING: Uncompress GTF annotation file or create from GFF3 if required
     */
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf, publish_genome ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff, publish_genome ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff, publish_genome ).gtf
        ch_software_versions = ch_software_versions.mix(GFFREAD.out.version)
    }

    /*
     * PREPROCESSING: Uncompress additional fasta file and concatenate with reference fasta and gtf files
     */
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( params.additional_fasta, publish_genome ).gunzip
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, publish_genome )
        ch_fasta = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf   = CAT_ADDITIONAL_FASTA.out.gtf
    }

    /*
     * PREPROCESSING: Uncompress BED12 annotation file or create from GTF if required
     */
    if (params.bed12) {
        if (params.bed12.endsWith('.gz')) {
            ch_bed12 = GUNZIP_BED12 ( params.bed12, publish_genome ).gunzip
        } else {
            ch_bed12 = file(params.bed12)
        }
    } else {
        ch_bed12 = GTF2BED ( ch_gtf, publish_genome )
    }

    /*
     * PREPROCESSING: Check genome/transcriptome indices and uncompress if required
     */
    def publish_index = params.save_reference ? [publish_dir : 'genome/index'] : [publish_files : [:]]
    if (!params.skip_alignment) {
        if (params.aligner == 'star') {
            if (params.star_index) {
                if (params.star_index.endsWith('.tar.gz')) {
                    ch_star_index = UNTAR_STAR_INDEX ( params.star_index, publish_index ).untar
                } else {
                    ch_star_index = file(params.star_index)
                }
            } else {
                // TODO nf-core: Not working - only save indices if --save_reference is specified
                if (params.save_reference) { params.modules['star_genomegenerate']['publish_files'] = null }
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf, params.modules['star_genomegenerate'] ).index
            }
        } else if (params.aligner == 'hisat2') {
            if (params.hisat2_index) {
                if (params.hisat2_index.endsWith('.tar.gz')) {
                    ch_hisat2_index = UNTAR_HISAT2_INDEX ( params.hisat2_index, publish_index ).untar
                } else {
                    ch_hisat2_index = file(params.hisat2_index)
                }
            } else {
                if (!params.splicesites) {
                    def publish_splicesites = params.save_reference ? [publish_dir : 'genome/index/hisat2'] : [publish_files : [:]]
                    ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf, publish_splicesites ).txt
                } else {
                    ch_splicesites = file(params.splicesites)
                }
                // TODO nf-core: Not working - only save indices if --save_reference is specified
                if (params.save_reference) { params.modules['hisat2_build']['publish_files'] = null }
                ch_hisat2_index = HISAT2_BUILD ( ch_fasta, ch_gtf, ch_splicesites, params.modules['hisat2_build'] ).index
            }
        }

        if (!skip_rsem && params.aligner == "star") {
            if (params.rsem_index) {
                if (params.rsem_index.endsWith('.tar.gz')) {
                    ch_rsem_index = UNTAR_RSEM_INDEX ( params.rsem_index, publish_index ).untar
                } else {
                    ch_rsem_index = file(params.rsem_index)
                }
            } else {
                // TODO nf-core: Not working - only save indices if --save_reference is specified
                if (params.save_reference) { params.modules['rsem_preparereference']['publish_files'] = null }
                //println(params.modules['rsem_preparereference'])
                ch_rsem_index = RSEM_PREPAREREFERENCE ( ch_fasta, ch_gtf, params.modules['rsem_preparereference'] ).index
            }
        }
    }

    if (params.pseudo_aligner == 'salmon') {
        if (params.salmon_index) {
            if (params.salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX ( params.salmon_index, publish_index ).untar
            } else {
                ch_salmon_index = file(params.salmon_index)
            }
        } else {
            if (params.transcript_fasta) {
                if (params.transcript_fasta.endsWith('.gz')) {
                    ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( params.transcript_fasta, publish_genome ).gunzip
                } else {
                    ch_transcript_fasta = file(params.transcript_fasta)
                }
            } else {
                def publish_transcripts = params.save_reference ? [publish_dir : 'genome/index/salmon'] : [publish_files : [:]]
                ch_transcript_fasta = TRANSCRIPTS2FASTA ( ch_fasta, ch_gtf, publish_transcripts ).fasta
            }
            // TODO nf-core: Not working - only save indices if --save_reference is specified
            if (params.save_reference) { params.modules['salmon_index']['publish_files'] = null }
            def gencode = params.gencode  ? "--gencode" : ""
            params.modules['salmon_index'].args += gencode
            ch_salmon_index = SALMON_INDEX ( ch_transcript_fasta, params.modules['salmon_index'])
        }
    }

    /*
     * Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input, params.seq_center, [:] )

    /*
     * Concatenate FastQ files from same sample if required
     */

    /*
     * Read QC & trimming
     */
    nextseq = params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
    params.modules['trimgalore'].args += nextseq
    if (params.save_trimmed) { params.modules['trimgalore'].publish_files.put('fq.gz','') }
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        params.skip_trimming,
        params.modules['fastqc'],
        params.modules['trimgalore']
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    ch_sortmerna_fasta = Channel.from(ch_ribo_db.readLines()).map { row -> file(row) }
//
//     /*
//      * Map reads & BAM QC
//      */
//     score = params.bwa_min_score ? " -T ${params.bwa_min_score}" : ''
//     params.modules['bwa_mem'].args += score
//     MAP_BWA_MEM (
//         FASTQC_TRIMGALORE.out.reads,
//         ch_index,
//         ch_fasta,
//         params.modules['bwa_mem'],
//         params.modules['samtools_sort_lib']
//     )
//     ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.bwa_version.first())
//     ch_software_versions = ch_software_versions.mix(MAP_BWA_MEM.out.samtools_version.first().ifEmpty(null))
//
//     /*
//      * Mark duplicates & filter BAM files
//      */
//     MARK_DUPLICATES_PICARD (
//         PICARD_MERGESAMFILES.out.bam,
//         params.modules['picard_markduplicates'],
//         params.modules['samtools_sort_merged_lib']
//     )
//
//     /*
//      * Post alignment QC
//      */
//     PICARD_COLLECTMULTIPLEMETRICS (
//         BAM_CLEAN.out.bam,
//         ch_fasta,
//         params.modules['picard_collectmultiplemetrics']
//     )
//
//     PRESEQ_LCEXTRAP (
//         BAM_CLEAN.out.bam,
//         params.modules['preseq_lcextrap']
//     )
//     ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))
//
//     /*
//      * Coverage tracks
//      */
//     BEDTOOLS_GENOMECOV (
//         BAM_CLEAN.out.bam.join(BAM_CLEAN.out.flagstat, by: [0]),
//         params.modules['bedtools_genomecov']
//     )
//
//     UCSC_BEDRAPHTOBIGWIG (
//         BEDTOOLS_GENOMECOV.out.bedgraph,
//         GET_CHROM_SIZES.out.sizes,
//         params.modules['ucsc_bedgraphtobigwig']
//     )
//     ch_software_versions = ch_software_versions.mix(UCSC_BEDRAPHTOBIGWIG.out.version.first().ifEmpty(null))
//
    /*
     * Pipeline reporting
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
//         PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]),
//
//         PRESEQ_LCEXTRAP.out.ccurve.collect{it[1]}.ifEmpty([]),
//
//         SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]),
//         // path ('macs/consensus/*') from ch_macs_consensus_deseq_mqc.collect().ifEmpty([])
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

// if (params.with_umi) {
//     process UMITOOLS_EXTRACT {
//         tag "$name"
//         label "low_memory"
//         publishDir "${params.outdir}/umitools/extract", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.endsWith('.log')) filename
//                 else if (!params.save_umi_intermeds && filename == "where_are_my_files.txt") filename
//                 else if (params.save_umi_intermeds && filename != "where_are_my_files.txt") filename
//                 else null
//             }
//
//         input:
//         tuple val(name), path(reads) from ch_raw_reads_umitools
//         path wherearemyfiles from ch_where_are_my_files
//
//         output:
//         tuple val(name), path("*fq.gz") into raw_reads_trimgalore
//         path "*.log"
//         path "where_are_my_files.txt"
//
//         script:
//         if (params.single_end) {
//             """
//             umi_tools extract \\
//                 -I $reads \\
//                 -S ${name}_umi_extracted.fq.gz \\
//                 --extract-method=${params.umitools_extract_method} \\
//                 --bc-pattern="${params.umitools_bc_pattern}" \\
//                 ${params.umitools_extract_extra} > ${name}_umi_extract.log
//             """
//         }  else {
//             """
//             umi_tools extract \\
//                 -I ${reads[0]} \\
//                 --read2-in=${reads[1]} \\
//                 -S ${name}_umi_extracted_R1.fq.gz \\
//                 --read2-out=${name}_umi_extracted_R2.fq.gz \\
//                 --extract-method=${params.umitools_extract_method} \\
//                 --bc-pattern="${params.umitools_bc_pattern}" \\
//                 ${params.umitools_extract_extra} > ${name}_umi_extract.log
//             """
//         }
//     }
// } else {
//     raw_reads_trimgalore = ch_raw_reads_umitools
//     umi_tools_extract_results = Channel.empty()
// }
//
//
// if (!params.remove_ribo_rna) {
//     trimgalore_reads
//         .into { trimmed_reads_alignment
//                 trimmed_reads_salmon }
//     sortmerna_logs = Channel.empty()
// } else {
//     process SORTMERNA {
//         tag "$name"
//         label 'low_memory'
//         publishDir "${params.outdir}/sortmerna", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("_rRNA_report.txt") > 0) "logs/$filename"
//                 else if (params.save_non_ribo_reads) "reads/$filename"
//                 else null
//             }
//
//         input:
//         tuple val(name), path(reads) from trimgalore_reads
//         path fasta from sortmerna_fasta.collect()
//
//         output:
//         tuple val(name), path("*.fq.gz") into trimmed_reads_alignment,
//                                               trimmed_reads_salmon
//         path "*_rRNA_report.txt" into sortmerna_logs
//
//         script:
//         //concatenate reference files: ${db_fasta},${db_name}:${db_fasta},${db_name}:...
//         def Refs = ""
//         for (i=0; i<fasta.size(); i++) { Refs+= " --ref ${fasta[i]}" }
//         if (params.single_end) {
//             """
//             sortmerna \\
//                 $Refs \\
//                 --reads $reads \\
//                 --num_alignments 1 \\
//                 --threads $task.cpus \\
//                 --workdir . \\
//                 --fastx \\
//                 --aligned rRNA-reads \\
//                 --other non-rRNA-reads \\
//                 -v
//
//             gzip --force < non-rRNA-reads.fq > ${name}.fq.gz
//             mv rRNA-reads.log ${name}_rRNA_report.txt
//             """
//         } else {
//             """
//             sortmerna \\
//                 $Refs \\
//                 --reads ${reads[0]} \\
//                 --reads ${reads[1]} \\
//                 --num_alignments 1 \\
//                 --threads $task.cpus \\
//                 --workdir . \\
//                 --fastx \\
//                 --aligned rRNA-reads \\
//                 --other non-rRNA-reads \\
//                 --paired_in \\
//                 --out2 \\
//                 -v
//
//             gzip --force < non-rRNA-reads_fwd.fq > ${name}_1.fq.gz
//             gzip --force < non-rRNA-reads_rev.fq > ${name}_2.fq.gz
//             mv rRNA-reads.log ${name}_rRNA_report.txt
//             """
//         }
//     }
// }
//
// // Function that checks the alignment rate of the STAR output
// // and returns true if the alignment passed and otherwise false
// good_alignment_scores = [:]
// poor_alignment_scores = [:]
// def check_log(logs) {
//     def percent_aligned = 0;
//     logs.eachLine { line ->
//         if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
//             percent_aligned = matcher[0][1]
//         }
//     }
//     logname = logs.getBaseName() - '.Log.final'
//     c_reset = params.monochrome_logs ? '' : "\033[0m";
//     c_green = params.monochrome_logs ? '' : "\033[0;32m";
//     c_red = params.monochrome_logs ? '' : "\033[0;31m";
//     if (percent_aligned.toFloat() <= params.percent_aln_skip.toFloat()) {
//         log.info "#${c_red}################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<${c_reset}"
//         poor_alignment_scores[logname] = percent_aligned
//         return false
//     } else {
//         log.info "-${c_green}           Passed alignment > star ($logname)   >> ${percent_aligned}% <<${c_reset}"
//         good_alignment_scores[logname] = percent_aligned
//         return true
//     }
// }
//
// if (!params.skip_alignment) {
//     if (params.aligner == 'star') {
//         hisat_stdout = Channel.empty()
//         process STAR_ALIGN {
//             tag "$name"
//             label 'high_memory'
//             publishDir "${params.outdir}/star", mode: params.publish_dir_mode,
//                 saveAs: { filename ->
//                     if (filename.indexOf(".bam") == -1) "logs/$filename"
//                     else if (params.save_unaligned && filename != "where_are_my_files.txt" && 'Unmapped' in filename) unmapped/filename
//                     else if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
//                     else if (params.save_align_intermeds && filename != "where_are_my_files.txt") filename
//                     else null
//                 }
//
//             input:
//             tuple val(name), path(reads) from trimmed_reads_alignment
//             path index from ch_star_index
//             path gtf from ch_gtf
//             path wherearemyfiles from ch_where_are_my_files
//
//             output:
//             tuple path("*Log.final.out"), path('*.sortedByCoord.out.bam'), path('*.toTranscriptome.out.bam') into star_aligned
//             path "*.out" into alignment_logs
//             path "*SJ.out.tab"
//             path "*Log.out" into star_log
//             path "where_are_my_files.txt"
//             path "*Unmapped*" optional true
//             path "${prefix}.Aligned.sortedByCoord.out.bam.bai" into bam_index
//
//             script:
//             prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
//             seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
//             def star_mem = task.memory ?: params.star_memory ?: false
//             def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
//             def unaligned = params.save_unaligned ? "--outReadsUnmapped Fastx" : ''
//             """
//             STAR \\
//                 --genomeDir $index \\
//                 --sjdbGTFfile $gtf \\
//                 --readFilesIn $reads  \\
//                 --runThreadN $task.cpus \\
//                 --twopassMode Basic \\
//                 --outWigType bedGraph \\
//                 --outSAMtype BAM SortedByCoordinate $avail_mem \\
//                 --readFilesCommand zcat \\
//                 --runDirPerm All_RWX $unaligned \\
//                 --quantMode TranscriptomeSAM \\
//                 --outFileNamePrefix $prefix. $seq_center \\
//                 --runRNGseed 0 \\
//                 $params.star_align_options
//
//             samtools index ${prefix}.Aligned.sortedByCoord.out.bam
//             """
//         }
//         // Filter removes all 'aligned' channels that fail the check
//         star_bams = Channel.create()
//         star_bams_transcriptome = Channel.create()
//         star_aligned
//             .filter { logs, bams, bams_transcriptome -> check_log(logs) }
//             .separate (star_bams, star_bams_transcriptome) {
//                 bam_set -> [bam_set[1], bam_set[2]]
//             }
//         bam = star_bams
//         bam_transcriptome = star_bams_transcriptome
//     }
//
//     if (params.aligner == 'hisat2') {
//         star_log = Channel.empty()
//         process HISAT2_ALIGN {
//             tag "$name"
//             label 'high_memory'
//             publishDir "${params.outdir}/hisat2", mode: params.publish_dir_mode,
//                 saveAs: { filename ->
//                     if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
//                     else if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
//                     else if (params.save_align_intermeds && filename != "where_are_my_files.txt") filename
//                     else null
//                 }
//
//             input:
//             tuple val(name), path(reads) from trimmed_reads_alignment
//             path index from ch_hisat2_index
//             path splicesites from ch_splicesites
//             path wherearemyfiles from ch_where_are_my_files
//
//             output:
//             path "${prefix}.bam" into hisat2_bam
//             path "${prefix}.hisat2_summary.txt" into alignment_logs
//             path "where_are_my_files.txt"
//             path "unmapped.hisat2*" optional true
//
//             script:
//             index_base = index[0].toString() - ~/.\d.ht2l?/
//             prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
//             seq_center = params.seq_center ? "--rg-id ${prefix} --rg CN:${params.seq_center.replaceAll('\\s','_')} SM:$prefix" : "--rg-id ${prefix} --rg SM:$prefix"
//             def rnastrandness = ''
//             if (forward_stranded && !unstranded) {
//                 rnastrandness = params.single_end ? '--rna-strandness F' : '--rna-strandness FR'
//             } else if (reverse_stranded && !unstranded) {
//                 rnastrandness = params.single_end ? '--rna-strandness R' : '--rna-strandness RF'
//             }
//             if (params.single_end) {
//                 unaligned = params.save_unaligned ? "--un-gz unmapped.hisat2.gz" : ''
//                 """
//                 hisat2 \\
//                     -x $index_base \\
//                     -U $reads \\
//                     $rnastrandness \\
//                     --known-splicesite-infile $splicesites \\
//                     -p $task.cpus $unaligned \\
//                     --met-stderr \\
//                     --new-summary \\
//                     --dta \\
//                     $params.hisat2_align_options \\
//                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
//                     | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
//                 """
//             } else {
//                 unaligned = params.save_unaligned ? "--un-conc-gz unmapped.hisat2.gz" : ''
//                 """
//                 hisat2 \\
//                     -x $index_base \\
//                     -1 ${reads[0]} \\
//                     -2 ${reads[1]} \\
//                     $rnastrandness \\
//                     --known-splicesite-infile $splicesites \\
//                     --no-mixed \\
//                     --no-discordant \\
//                     -p $task.cpus $unaligned \\
//                     --met-stderr \\
//                     --new-summary \\
//                     $params.hisat2_align_options \\
//                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
//                     | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
//                 """
//             }
//         }
//
//         process HISAT2_SORT_BAM {
//             tag "${bam.baseName}"
//             label 'mid_memory'
//             publishDir "${params.outdir}/hisat2", mode: params.publish_dir_mode,
//                 saveAs: { filename ->
//                     if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
//                     else if (params.save_align_intermeds && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
//                     else null
//                 }
//
//             input:
//             path bam from hisat2_bam
//             path wherearemyfiles from ch_where_are_my_files
//
//             output:
//             path "${bam.baseName}.sorted.bam" into bam
//             path "${bam.baseName}.sorted.bam.bai" into bam_index
//             path "where_are_my_files.txt"
//
//             script:
//             def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
//             def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
//             """
//             samtools sort \\
//                 $bam \\
//                 -@ $task.cpus $avail_mem \\
//                 -o ${bam.baseName}.sorted.bam
//             samtools index ${bam.baseName}.sorted.bam
//             """
//         }
//     }
//
//     if (params.with_umi) {
//         // preseq does not work on deduplicated BAM file. Pass it the raw BAM file.
//         bam
//             .into { bam_umitools_dedup
//                     bam_preseq }
//         bam_index_umitools_dedup = bam_index
//
//         process UMITOOLS_DEDUP {
//             tag "${bam.baseName}"
//             label "mid_memory"
//             publishDir "${params.outdir}/umitools/dedup", mode: params.publish_dir_mode,
//                 saveAs: { filename ->
//                     if (filename.endsWith('.tsv')) filename
//                     else if (!params.save_umi_intermeds && filename == "where_are_my_files.txt") filename
//                     else if (params.save_umi_intermeds && filename != "where_are_my_files.txt") filename
//                     else null
//                 }
//
//             input:
//             path bam from bam_umitools_dedup
//             path bai from bam_index_umitools_dedup
//             path wherearemyfiles from ch_where_are_my_files
//
//             output:
//             path "*.bam" into bam_dedup
//             path "*.bai" into bam_dedup_index
//             path "where_are_my_files.txt"
//             path "*.tsv"
//
//             script:
//             """
//             umi_tools dedup \\
//                 -I $bam \\
//                 -S ${bam.baseName}_deduplicated.bam \\
//                 --output-stats=${bam.baseName} \\
//                 $params.umitools_dedup_extra
//             samtools index ${bam.baseName}_deduplicated.bam
//             """
//         }
//
//         // RSEM transcriptome BAM file treated separately...
//         if (!skip_rsem) {
//             process UMITOOLS_DEDUP_TRANSCRIPTOME {
//                 tag "${bam.baseName}"
//                 label "mid_memory"
//                 publishDir "${params.outdir}/umitools/dedup/transcriptome", mode: params.publish_dir_mode,
//                     saveAs: { filename ->
//                         if (filename.endsWith('.tsv')) filename
//                         else if (params.save_umi_intermeds) filename
//                         else null
//                     }
//
//                 input:
//                 path bam from bam_transcriptome
//
//                 output:
//                 path "*_deduplicated.bam" into bam_rsem
//                 path "*.tsv"
//
//                 script:
//                 // the transcriptome BAM file is not sorted or indexed by STAR
//                 // since this is the only process consuming this BAM file,
//                 // sorting and indexing happens right here.
//                 def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
//                 def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
//                 """
//                 samtools sort \\
//                     $bam \\
//                     -@ $task.cpus $avail_mem \\
//                     -o ${bam.baseName}.sorted.bam
//                 samtools index ${bam.baseName}.sorted.bam
//
//                 umi_tools dedup \\
//                     -I ${bam.baseName}.sorted.bam \\
//                     -S ${bam.baseName}_deduplicated.bam \\
//                     --output-stats=${bam.baseName} \\
//                     $params.umitools_dedup_extra
//                 """
//             }
//         }
//
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
//     process PRESEQ {
//         tag "${bam.baseName - '.sorted'}"
//         label 'high_time'
//         publishDir "${params.outdir}/preseq", mode: params.publish_dir_mode
//
//         when:
//         !params.skip_qc && !params.skip_preseq
//
//         input:
//         path bam from bam_preseq
//
//         output:
//         path "${bam.baseName}.ccurve.txt" into preseq_results
//
//         script:
//         """
//         preseq lc_extrap -v -B $bam -o ${bam.baseName}.ccurve.txt
//         """
//     }
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
// if (params.pseudo_aligner == 'salmon') {
//     process SALMON_QUANT {
//         tag "$sample"
//         publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
//
//         input:
//         tuple val(sample), path(reads) from trimmed_reads_salmon
//         path index from ch_salmon_index
//         path gtf from ch_gtf
//
//         output:
//         tuple val(sample), path("${sample}/") into salmon_parsegtf,
//                                                    salmon_tximport
//         path "${sample}/" into salmon_logs
//
//         script:
//         def rnastrandness = params.single_end ? 'U' : 'IU'
//         if (forward_stranded && !unstranded) {
//             rnastrandness = params.single_end ? 'SF' : 'ISF'
//         } else if (reverse_stranded && !unstranded) {
//             rnastrandness = params.single_end ? 'SR' : 'ISR'
//         }
//         def endedness = params.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
//         def unmapped = params.save_unaligned ? "--writeUnmappedNames" : ''
//         """
//         salmon quant \\
//             --validateMappings \\
//             --seqBias \\
//             --useVBOpt \\
//             --gcBias \\
//             --geneMap $gtf \\
//             --threads $task.cpus \\
//             --libType=$rnastrandness \\
//             --index $index \\
//             $endedness \\
//             $unmapped \\
//             -o $sample
//         """
//     }
//
//     process SALMON_TX2GENE {
//         label 'low_memory'
//         publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
//
//         input:
//         path ("salmon/*") from salmon_parsegtf.collect{it[1]}
//         path gtf from ch_gtf
//
//         output:
//         path "tx2gene.csv" into salmon_tx2gene,
//                                 salmon_merge_tx2gene
//
//         script:
//         """
//         parse_gtf.py --gtf $gtf --salmon salmon --id $params.fc_group_features --extra $params.fc_extra_attributes -o tx2gene.csv
//         """
//     }
//
//     process SALMON_TXIMPORT {
//         label 'low_memory'
//         publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
//
//         input:
//         tuple val(name), path("salmon/*") from salmon_tximport
//         path tx2gene from salmon_tx2gene.collect()
//
//         output:
//         path "${name}_salmon_gene_tpm.csv" into salmon_gene_tpm
//         path "${name}_salmon_gene_counts.csv" into salmon_gene_counts
//         path "${name}_salmon_transcript_tpm.csv" into salmon_transcript_tpm
//         path "${name}_salmon_transcript_counts.csv" into salmon_transcript_counts
//
//         script:
//         """
//         tximport.r NULL salmon $name
//         """
//     }
//
//     process SALMON_MERGE {
//         label 'mid_memory'
//         publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
//
//         input:
//         path gene_tpm_files from salmon_gene_tpm.collect()
//         path gene_count_files from salmon_gene_counts.collect()
//         path transcript_tpm_files from salmon_transcript_tpm.collect()
//         path transcript_count_files from salmon_transcript_counts.collect()
//         path tx2gene from salmon_merge_tx2gene
//
//         output:
//         path "salmon_merged*.csv" into salmon_merged_ch
//         path "*.rds"
//
//         script:
//         // First field is the gene/transcript ID
//         gene_ids = "<(cut -f1 -d, ${gene_tpm_files[0]} | tail -n +2 | cat <(echo '${params.fc_group_features}') - )"
//         transcript_ids = "<(cut -f1 -d, ${transcript_tpm_files[0]} | tail -n +2 | cat <(echo 'transcript_id') - )"
//
//         // Second field is counts/TPM
//         gene_tpm = gene_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
//         gene_counts = gene_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
//         transcript_tpm = transcript_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
//         transcript_counts = transcript_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
//         """
//         paste -d, $gene_ids $gene_tpm > salmon_merged_gene_tpm.csv
//         paste -d, $gene_ids $gene_counts > salmon_merged_gene_counts.csv
//         paste -d, $transcript_ids $transcript_tpm > salmon_merged_transcript_tpm.csv
//         paste -d, $transcript_ids $transcript_counts > salmon_merged_transcript_counts.csv
//
//         se.r NULL salmon_merged_gene_counts.csv salmon_merged_gene_tpm.csv
//         se.r NULL salmon_merged_transcript_counts.csv salmon_merged_transcript_tpm.csv
//         """
//     }
// } else {
//     salmon_logs = Channel.empty()
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
