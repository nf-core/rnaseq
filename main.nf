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
params.gene_bed     = params.genomes[ params.genome ]?.bed12
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
if (params.gtf && params.gff)   { Checks.gtf_gff_warn(log) }

// Check input path parameters to see if they exist
checkPathParamList = [
    params.gtf, params.gff, params.gene_bed, params.additional_fasta, params.transcript_fasta,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index,
    params.splicesites, params.ribo_database_manifest
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check rRNA databases for sortmerna
ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
if (ch_ribo_db.isEmpty()) {exit 1, "File ${ch_ribo_db.getName()} is empty!"}

// Check alignment parameters
def alignerList       = ['star', 'hisat2', 'star_rsem']
def pseudoAlignerList = ['salmon']
if (!params.skip_alignment) {
    if (!alignerList.contains(params.aligner)) {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: ${alignerList.join(', ')}"
    }
} else {
    if (!params.pseudo_aligner) {
        exit 1, "--skip_alignment specified without --pseudo_aligner...please specify e.g. --pseudo_aligner ${pseudoAlignerList[0]}"
    }
    Checks.skip_alignment_warn(log)
}
if (params.pseudo_aligner) {
    if (!pseudoAlignerList.contains(params.pseudo_aligner)) {
        exit 1, "Invalid pseudo aligner option: ${params.pseudo_aligner}. Valid options: ${pseudoAlignerList.join(', ')}"
    } else {
        if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
            exit 1, "To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff"
        }
    }
}

// Check and exit if we are using '--aligner star_rsem' and '--with_umi'
if (!params.skip_alignment && params.aligner == 'star_rsem' && params.with_umi) {
    Checks.rsem_umi_error(log)
    exit 1
}

// Check which RSeQC modules we are running
def rseqcModuleList = [
    'bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation',
    'junction_saturation', 'read_distribution', 'read_duplication'
]
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if ((rseqcModuleList + rseqc_modules).unique().size() != rseqcModuleList.size()) {
    exit 1, "Invalid RSeqC module options: ${params.rseqc_modules}. Valid options: ${rseqcModuleList.join(', ')}"
}

// Show a big warning message if we are using GRCh38 NCBI assembly
def skip_biotype_qc = params.skip_biotype_qc
if (params.genome == 'GRCh38') {
    if (params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
        Checks.genome_warn(log)
        skip_biotype_qc = true
    }
}

// Save AWS IGenomes file containing annotation version
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

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

// Header files for MultiQC
ch_mdsplot_header        = file("$baseDir/assets/multiqc/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header        = file("$baseDir/assets/multiqc/heatmap_header.txt", checkIfExists: true)
ch_biotypes_header       = file("$baseDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

summary = Schema.params_summary(workflow, params)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { CAT_FASTQ                   } from './modules/local/process/cat_fastq'
include { MULTIQC_CUSTOM_BIOTYPE      } from './modules/local/process/multiqc_custom_biotype'
include { MULTIQC_CUSTOM_FAIL_MAPPED  } from './modules/local/process/multiqc_custom_fail_mapped'
include { MULTIQC_CUSTOM_STRAND_CHECK } from './modules/local/process/multiqc_custom_strand_check'
include { FEATURECOUNTS_MERGE_COUNTS  } from './modules/local/process/featurecounts_merge_counts'
include { EDGER_CORRELATION           } from './modules/local/process/edger_correlation'
include { DUPRADAR                    } from './modules/local/process/dupradar'
include { GET_SOFTWARE_VERSIONS       } from './modules/local/process/get_software_versions'
include { MULTIQC                     } from './modules/local/process/multiqc'

include { INPUT_CHECK                 } from './modules/local/subworkflow/input_check'
include { PREPARE_GENOME              } from './modules/local/subworkflow/prepare_genome'
include { ALIGN_STAR                  } from './modules/local/subworkflow/align_star'
include { ALIGN_HISAT2                } from './modules/local/subworkflow/align_hisat2'
include { QUANTIFY_RSEM               } from './modules/local/subworkflow/quantify_rsem'
include { QUANTIFY_SALMON             } from './modules/local/subworkflow/quantify_salmon'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { SAMTOOLS_INDEX              } from './modules/nf-core/software/samtools/index/main'
include { PRESEQ_LCEXTRAP             } from './modules/nf-core/software/preseq/lcextrap/main'
include { SORTMERNA                   } from './modules/nf-core/software/sortmerna/main'
include { STRINGTIE                   } from './modules/nf-core/software/stringtie/main'
include { QUALIMAP_RNASEQ             } from './modules/nf-core/software/qualimap/rnaseq/main'
include { UMITOOLS_DEDUP              } from './modules/nf-core/software/umitools/dedup/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS
          SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_BIOTYPE } from './modules/nf-core/software/subread/featurecounts/main'

include { FASTQC_UMITOOLS_TRIMGALORE  } from './modules/nf-core/subworkflow/fastqc_umitools_trimgalore'
include { BAM_SORT_SAMTOOLS           } from './modules/nf-core/subworkflow/bam_sort_samtools'
include { MARK_DUPLICATES_PICARD      } from './modules/nf-core/subworkflow/mark_duplicates_picard'
include { RSEQC                       } from './modules/nf-core/subworkflow/rseqc'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []
def pass_percent_mapped = [:]
def fail_percent_mapped = [:]

workflow {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    def publish_genome_options  = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
    def publish_index_options   = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]
    if (!params.save_reference) { params.modules['gffread']['publish_files'] = false }
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.gene_bed,
        params.additional_fasta,
        params.modules['gffread'],
        publish_genome_options
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input, [:] )
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
    if (!params.save_merged_fastq) { params.modules['cat_fastq'].publish_files = false }
    CAT_FASTQ ( ch_cat_fastq, params.modules['cat_fastq'] )

    /*
     * SUBWORKFLOW: Read QC, extract UMI and trim adapters
     */
    def method  = params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : ''
    def pattern = params.umitools_bc_pattern     ? "--bc-pattern='${params.umitools_bc_pattern}'"       : ''
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
    ch_genome_bam        = Channel.empty()
    ch_genome_bai        = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_star_multiqc      = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star') {
        if (!params.save_reference)      { params.modules['star_genomegenerate']['publish_files'] = false }
        if (params.save_unaligned)       { params.modules['star_align'].publish_files.put('fastq.gz','unmapped') }
        def unaligned = params.save_unaligned ? " --outReadsUnmapped Fastx" : ''
        params.modules['star_align'].args += unaligned
        if (params.save_align_intermeds) { params.modules['star_align'].publish_files.put('bam','') }
        if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
            params.modules['samtools_sort'].publish_files.put('bam','')
            params.modules['samtools_sort'].publish_files.put('bai','')
        }
        
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
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
     */
    ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        if (!params.save_reference) { params.modules['rsem_preparereference']['publish_files'] = false }
        if (params.save_align_intermeds) { params.modules['rsem_calculateexpression'].publish_files.put('bam','') }
        if (params.save_align_intermeds || params.skip_markduplicates) {
            params.modules['samtools_sort'].publish_files.put('bam','')
            params.modules['samtools_sort'].publish_files.put('bai','')
        }

        QUANTIFY_RSEM (
            ch_trimmed_reads,
            params.rsem_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            publish_index_options,
            params.modules['rsem_preparereference'],
            params.modules['rsem_calculateexpression'],
            params.modules['samtools_sort'],
            params.modules['rsem_merge_counts']
        )
        ch_genome_bam        = QUANTIFY_RSEM.out.bam
        ch_genome_bai        = QUANTIFY_RSEM.out.bai
        ch_samtools_stats    = QUANTIFY_RSEM.out.stats
        ch_samtools_flagstat = QUANTIFY_RSEM.out.flagstat
        ch_samtools_idxstats = QUANTIFY_RSEM.out.idxstats
        ch_star_multiqc      = QUANTIFY_RSEM.out.logs
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
        ch_software_versions = ch_software_versions.mix(QUANTIFY_RSEM.out.rsem_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_RSEM.out.samtools_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Alignment with HISAT2
     */
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        if (!params.save_reference)      { params.modules['hisat2_build']['publish_files'] = false }
        if (params.save_unaligned)       { params.modules['hisat2_align'].publish_files.put('fastq.gz','unmapped') }
        if (params.save_align_intermeds) { params.modules['hisat2_align'].publish_files.put('bam','') }
        if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
            params.modules['samtools_sort'].publish_files.put('bam','')
            params.modules['samtools_sort'].publish_files.put('bai','')
        }

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
     * Filter channels to get samples that passed STAR minimum mapping percentage
     */
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + Checks.get_star_percent_mapped(workflow, params, log, align_log) }
            .set { ch_percent_mapped }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bai
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bai }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_fail_mapping_multiqc = MULTIQC_CUSTOM_FAIL_MAPPED ( ch_pass_fail_mapped.fail.collect(), [publish_files: false] )
    }

    /*
     * MODULE: Run Preseq
     */
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP ( ch_genome_bam, params.modules['preseq_lcextrap'] )
        ch_preseq_multiqc    = PRESEQ_LCEXTRAP.out.ccurve
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.version.first().ifEmpty(null))
    }

    /*
     * MODULE: Remove duplicate reads from BAM file based on UMIs
     */
    if (!params.skip_alignment && params.aligner != 'star_rsem' && params.with_umi) {
        if (params.save_align_intermeds || params.skip_markduplicates || params.save_umi_intermeds) {
            params.modules['umitools_dedup'].publish_files.put('bam','')
            params.modules['umitools_dedup'].publish_files.put('bai','')
        }

        UMITOOLS_DEDUP ( ch_genome_bam.join(ch_genome_bai, by: [0]), params.modules['umitools_dedup'] )
        SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam, params.modules['umitools_dedup'] )
        ch_genome_bam = UMITOOLS_DEDUP.out.bam
        ch_genome_bai = SAMTOOLS_INDEX.out.bai
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
     * MODULE: STRINGTIE
     */
    if (!params.skip_alignment && !params.skip_stringtie) {
        def ignore_gtf = params.stringtie_ignore_gtf ? "" : " -e"
        params.modules['stringtie'].args += ignore_gtf
        STRINGTIE ( ch_genome_bam, PREPARE_GENOME.out.gtf, params.modules['stringtie'] )
        ch_software_versions = ch_software_versions.mix(STRINGTIE.out.version.first().ifEmpty(null))
    }

    /*
     * MODULE: Count reads relative to features using featureCounts
     */
    ch_edger_multiqc = Channel.empty()
    ch_featurecounts_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner != 'star_rsem' && !params.skip_featurecounts) {
        def fc_extra_attributes = params.fc_extra_attributes  ? " --extraAttributes $params.fc_extra_attributes" : ""
        params.modules['subread_featurecounts'].args += fc_extra_attributes
        params.modules['subread_featurecounts'].args += " -g $params.fc_group_features -t $params.fc_count_type"

        SUBREAD_FEATURECOUNTS ( ch_genome_bam.combine(PREPARE_GENOME.out.gtf), params.modules['subread_featurecounts'] )
        ch_featurecounts_multiqc = SUBREAD_FEATURECOUNTS.out.summary
        ch_software_versions = ch_software_versions.mix(SUBREAD_FEATURECOUNTS.out.version.first().ifEmpty(null))

        FEATURECOUNTS_MERGE_COUNTS ( SUBREAD_FEATURECOUNTS.out.counts.collect{it[1]}, params.modules['featurecounts_merge_counts'] )

        if (!params.skip_qc & !params.skip_edger) {    
            EDGER_CORRELATION ( SUBREAD_FEATURECOUNTS.out.counts.collect{it[1]}, ch_mdsplot_header, ch_heatmap_header, params.modules['edger_correlation'] )
            ch_edger_multiqc = EDGER_CORRELATION.out.multiqc
            ch_software_versions = ch_software_versions.mix(EDGER_CORRELATION.out.version.ifEmpty(null))
        }
    }

    /*
     * MODULE: Feature biotype QC using featureCounts
     */
    ch_featurecounts_biotype_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_featurecounts && !skip_biotype_qc) {
        def biotype = params.gencode ? "gene_type" : params.fc_group_features_type
        params.modules['subread_featurecounts_biotype'].args += " -g $biotype -t $params.fc_count_type"
        SUBREAD_FEATURECOUNTS_BIOTYPE ( ch_genome_bam.combine(PREPARE_GENOME.out.gtf), params.modules['subread_featurecounts_biotype'] )

        MULTIQC_CUSTOM_BIOTYPE ( SUBREAD_FEATURECOUNTS_BIOTYPE.out.counts, ch_biotypes_header, params.modules['multiqc_custom_biotype'] )
        ch_featurecounts_biotype_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
    }

    /*
     * MODULE: Downstream QC steps
     */
    ch_qualimap_multiqc           = Channel.empty()
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {
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
        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            RSEQC (
                ch_genome_bam,
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules,
                params.modules['rseqc_bamstat'],
                params.modules['rseqc_innerdistance'],
                params.modules['rseqc_inferexperiment'],
                params.modules['rseqc_junctionannotation'],
                params.modules['rseqc_junctionsaturation'],
                params.modules['rseqc_readdistribution'],
                params.modules['rseqc_readduplication']
            )
            ch_bamstat_multiqc            = RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = RSEQC.out.readduplication_pos_xls
            ch_software_versions = ch_software_versions.mix(RSEQC.out.version.first().ifEmpty(null))

            ch_inferexperiment_multiqc
                .map { meta, strand_log -> [ meta ] + Checks.get_inferexperiment_strandedness(strand_log, 30) }
                .filter { it[0].strandedness != it[1] }
                .map { meta, strandedness, sense, antisense, undetermined ->
                    [ "$meta.id\t$meta.strandedness\t$strandedness\t$sense\t$antisense\t$undetermined" ]
                }
                .set { ch_fail_strand }

            ch_fail_strand_multiqc = MULTIQC_CUSTOM_STRAND_CHECK ( ch_fail_strand.collect(), [publish_files: false] )
        }
    }

    /*
     * SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
     */
    ch_salmon_multiqc = Channel.empty()
    if (params.pseudo_aligner == 'salmon') {
        if (!params.save_reference) { params.modules['salmon_index']['publish_files'] = false }
        def gencode = params.gencode  ? " --gencode" : ""
        params.modules['salmon_index'].args += gencode

        def unmapped = params.save_unaligned ? " --writeUnmappedNames" : ''
        params.modules['salmon_quant'].args += unmapped

        QUANTIFY_SALMON (
            ch_trimmed_reads,
            params.salmon_index,
            params.transcript_fasta,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            publish_index_options,
            publish_genome_options,
            params.modules['salmon_index'],
            params.modules['salmon_quant'],
            params.modules['salmon_merge_counts']
        )
        ch_salmon_multiqc    = QUANTIFY_SALMON.out.results
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.salmon_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.tximeta_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.summarizedexperiment_version.ifEmpty(null))
    }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect(), [publish_files : ['csv':'']] )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary     = Schema.params_summary_multiqc(summary)
        ch_workflow_summary  = Channel.value(workflow_summary)

        if (params.skip_alignment) { params.modules['multiqc']['publish_dir'] = '' }
        def multiqc_title = params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''
        params.modules['multiqc'].args += "$multiqc_title"
        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_mapping_multiqc.ifEmpty([]),
            ch_fail_strand_multiqc.ifEmpty([]),
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
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_edger_multiqc.collect().ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_biotype_multiqc.collect{it[1]}.ifEmpty([]),
            params.modules['multiqc']
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, summary, baseDir, multiqc_report, log, fail_percent_mapped)
    Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
