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

// Auto-load genome files from genome config
params.fasta        = params.genomes[ params.genome ]?.fasta
params.gtf          = params.genomes[ params.genome ]?.gtf
params.gff          = params.genomes[ params.genome ]?.gff
params.bed12        = params.genomes[ params.genome ]?.bed12
params.star_index   = params.genomes[ params.genome ]?.star
params.hisat2_index = params.genomes[ params.genome ]?.hisat2
params.rsem_index   = params.genomes[ params.genome ]?.rsem
params.salmon_index = params.genomes[ params.genome ]?.salmon
anno_readme         = params.genomes[ params.genome ]?.readme

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check mandatory parameters
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, 'Genome fasta file not specified!' }
if (!params.gtf && !params.gff) { exit 1, "No GTF or GFF3 annotation specified!" }
if (params.gtf && params.gff)   { log.info "WARN: Both GTF and GFF have been provided: Using GTF as priority." }

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
if (!params.skip_alignment) {
    if (params.aligner != 'star' && params.aligner != 'hisat2') {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
    }
    if (params.aligner != "star" && !params.skip_rsem) {
        log.info "WARN: RSEM only works when '--aligner star' is set. Disabling RSEM."
    }
} else {
    log.info "WARN: Skipping alignment processes..."
    if (!params.pseudo_aligner) {
        exit 1, "--skip_alignment specified without --pseudo_aligner .. did you mean to specify --pseudo_aligner salmon"
    }
}
if (params.pseudo_aligner) {
    if (params.pseudo_aligner != 'salmon') {
        exit 1, "Invalid pseudo aligner option: ${params.pseudo_aligner}. Valid options: 'salmon'"
    } //else {
        // if (!params.salmon_index || !params.transcript_fasta || (!params.fasta || !(params.gff || params.gtf))) {
        //     exit 1, "To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff"
        // }
    //}
}

// Save AWS IGenomes file containing annotation version
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
 * Initiate parameters
 */
// Set biotype for featureCounts
def biotype = params.gencode ? "gene_type" : params.fc_group_features_type

/*
 * Check other parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

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

include { CAT_FASTQ                                      } from './modules/local/process/cat_fastq'
include { SORTMERNA                                      } from './modules/local/process/sortmerna'
include { UMITOOLS_DEDUP as UMITOOLS_DEDUP_GENOME
          UMITOOLS_DEDUP as UMITOOLS_DEDUP_TRANSCRIPTOME } from './modules/local/process/umitools_dedup'
// include { STRINGTIE                                      } from './modules/local/process/stringtie'
include { FEATURECOUNTS_MERGE_COUNTS                     } from './modules/local/process/featurecounts_merge_counts'
// include { EDGER_CORRELATION                             } from './modules/local/process/edger_correlation'
// include { RSEQC                                          } from './modules/local/process/rseqc'
include { QUALIMAP_RNASEQ                                } from './modules/local/process/qualimap_rnaseq'
include { DUPRADAR                                       } from './modules/local/process/dupradar'
include { OUTPUT_DOCUMENTATION                           } from './modules/local/process/output_documentation'
include { GET_SOFTWARE_VERSIONS                          } from './modules/local/process/get_software_versions'
include { MULTIQC                                        } from './modules/local/process/multiqc'

include { INPUT_CHECK                                    } from './modules/local/subworkflow/input_check'
include { PREPARE_GENOME                                 } from './modules/local/subworkflow/prepare_genome'
include { ALIGN_STAR                                     } from './modules/local/subworkflow/align_star'
include { ALIGN_HISAT2                                   } from './modules/local/subworkflow/align_hisat2'
include { QUANTIFY_RSEM                                  } from './modules/local/subworkflow/quantify_rsem'
include { QUANTIFY_SALMON                                } from './modules/local/subworkflow/quantify_salmon'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { SAMTOOLS_INDEX             } from './modules/nf-core/software/samtools/index/main'
include { PRESEQ_LCEXTRAP            } from './modules/nf-core/software/preseq/lcextrap/main'
include { SUBREAD_FEATURECOUNTS      } from './modules/nf-core/software/subread/featurecounts/main'

include { FASTQC_UMITOOLS_TRIMGALORE } from './modules/nf-core/subworkflow/fastqc_umitools_trimgalore'
include { BAM_SORT_SAMTOOLS          } from './modules/nf-core/subworkflow/bam_sort_samtools'
include { MARK_DUPLICATES_PICARD     } from './modules/nf-core/subworkflow/mark_duplicates_picard'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    def publish_genome_options = params.save_reference ? [publish_dir : 'genome']       : [publish_files : [:]]
    def publish_index_options  = params.save_reference ? [publish_dir : 'genome/index'] : [publish_files : [:]]
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.bed12,
        params.additional_fasta,
        publish_genome_options
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

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
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        if (params.save_non_ribo_reads) { params.modules['sortmerna'].publish_files.put('fastq.gz','') }
        ch_sortmerna_fasta = Channel.from(ch_ribo_db.readLines()).map { row -> file(row) }.collect()
        SORTMERNA ( ch_trimmed_reads, ch_sortmerna_fasta, params.modules['sortmerna'] )
            .reads
            .set { ch_trimmed_reads }
        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_software_versions = ch_software_versions.mix(SORTMERNA.out.version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Alignment with STAR
     */
    def good_alignment_scores = [:]
    def poor_alignment_scores = [:]
    ch_genome_bam        = Channel.empty()
    ch_genome_bai        = Channel.empty()
    ch_transcriptome_bam = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_star_multiqc      = Channel.empty()
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
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            params.modules['star_genomegenerate'],
            params.modules['star_align'],
            params.modules['samtools_sort']
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bai        = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))

        // Filter channels to get samples that passed minimum mapping percentage
        ch_star_multiqc
            .map { meta, align_log ->
                def percent_aligned = Checks.get_star_percent_mapped(params, log, align_log)
                if (percent_aligned <= params.percent_aln_skip.toFloat()) {
                    poor_alignment_scores[meta.id] = percent_aligned
                    [ meta, true ]
                } else {
                    good_alignment_scores[meta.id] = percent_aligned
                    [ meta, false ]
                }
            }
            .set { ch_sample_filter }

        ch_genome_bam
            .join(ch_sample_filter, by: [0])
            .map { meta, ofile, filter -> if (!filter) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bai
            .join(ch_sample_filter, by: [0])
            .map { meta, ofile, filter -> if (!filter) [ meta, ofile ] }
            .set { ch_genome_bai }

        ch_transcriptome_bam
            .join(ch_sample_filter, by: [0])
            .map { meta, ofile, filter -> if (!filter) [ meta, ofile ] }
            .set { ch_transcriptome_bam }
    }

    /*
     * SUBWORKFLOW: Alignment with HISAT2
     */
    ch_hisat2_multiqc = Channel.empty()
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
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            params.splicesites,
            params.modules['hisat2_build'],
            params.modules['hisat2_align'],
            params.modules['samtools_sort']
        )
        ch_genome_bam        = ALIGN_HISAT2.out.bam
        ch_genome_bai        = ALIGN_HISAT2.out.bai
        ch_samtools_stats    = ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = ALIGN_HISAT2.out.summary
        ch_software_versions = ch_software_versions.mix(ALIGN_HISAT2.out.hisat2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_HISAT2.out.samtools_version.first().ifEmpty(null))
    }

    /*
     * MODULE: Run Preseq
     */
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP ( ch_genome_bam, params.modules['preseq_lcextrap'] )
        ch_preseq_multiqc    = PRESEQ_LCEXTRAP.out.ccurve
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))
    }

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

        if (params.aligner == 'star' && !params.skip_rsem) {
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
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_genome_bam,
            params.modules['picard_markduplicates'],
            params.modules['picard_markduplicates_samtools']
        )
        ch_genome_bam             = MARK_DUPLICATES_PICARD.out.bam
        ch_genome_bai             = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Gene/transcript quantification with RSEM
     */
    ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star' && !params.skip_rsem) {
        // TODO nf-core: Not working - only save indices if --save_reference is specified
        if (params.save_reference) { params.modules['rsem_preparereference']['publish_files'] = null }
        QUANTIFY_RSEM (
            ch_transcriptome_bam,
            params.rsem_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            params.modules['rsem_preparereference'],
            params.modules['rsem_calculateexpression'],
            params.modules['rsem_merge_counts']
        )
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
        ch_software_versions = ch_software_versions.mix(QUANTIFY_RSEM.out.version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
     */
    ch_salmon_multiqc = Channel.empty()
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
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            publish_index_options,
            publish_genome_options,
            params.modules['salmon_index'],
            params.modules['salmon_quant'],
            params.modules['salmon_merge_counts']
        )
        ch_salmon_multiqc    = QUANTIFY_SALMON.out.results
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.version.first().ifEmpty(null))
    }

    /*
     * MODULE: Downstream QC steps
     */
    ch_rseqc_multiqc    = Channel.empty()
    ch_qualimap_multiqc = Channel.empty()
    ch_dupradar_multiqc = Channel.empty()
    ch_edger_multiqc    = Channel.empty()
    if (!params.skip_qc) {
        // if (!params.skip_rseqc) {
        //     RSEQC ( ch_genome_bam.join(ch_genome.bai, by: [0]), PREPARE_GENOME.out.bed12, params.modules['rseqc'] )
        //     ch_software_versions = ch_software_versions.mix(RSEQC.out.version.first().ifEmpty(null))
        // }
        if (!params.skip_qualimap) {
            QUALIMAP_RNASEQ ( ch_genome_bam, PREPARE_GENOME.out.gtf, params.modules['qualimap_rnaseq'] )
            ch_qualimap_multiqc  = QUALIMAP_RNASEQ.out.results
            ch_software_versions = ch_software_versions.mix(QUALIMAP_RNASEQ.out.version.first().ifEmpty(null))
        }
        if (!params.skip_dupradar) {
            DUPRADAR ( ch_genome_bam, PREPARE_GENOME.out.gtf, params.modules['dupradar'] )
            ch_dupradar_multiqc  = DUPRADAR.out.multiqc
            ch_software_versions = ch_software_versions.mix(DUPRADAR.out.version.first().ifEmpty(null))
        }
        // if (!params.skip_edger) {
        //     EDGER_CORRELATION ( counts.collect{it[1]}, ch_mdsplot_header, ch_heatmap_header, params.modules['edger_correlation'] )
        //     ch_edger_multiqc     = EDGER_CORRELATION.out.multiqc
        //     ch_software_versions = ch_software_versions.mix(EDGER_CORRELATION.out.version.first().ifEmpty(null))
        // }
    }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect(), [publish_files : ['csv':'']] )
    OUTPUT_DOCUMENTATION  ( ch_output_docs, ch_output_docs_images, [:] )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_mqc_summary(summary)
        ch_workflow_summary = Channel.value(workflow_summary)
        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rseqc_multiqc.collect{it[1]}.ifEmpty([]),                 // rseqc_results.collect().ifEmpty([])
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            // SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]) // featureCounts_logs.collect().ifEmpty([])
            // path ('featurecounts/biotype/*')                             // featureCounts_biotype.collect().ifEmpty([])
            ch_edger_multiqc.collect().ifEmpty([]),                      // sample_correlation_results.collect().ifEmpty([])
            params.modules['multiqc']
        )
    }
}

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
