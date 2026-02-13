#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.gene_bed         = getGenomeAttribute('bed12')
params.bbsplit_index    = getGenomeAttribute('bbsplit')
params.sortmerna_index  = getGenomeAttribute('sortmerna')
params.star_index       = getGenomeAttribute('star')
params.rsem_index       = getGenomeAttribute('rsem')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.salmon_index     = getGenomeAttribute('salmon')
params.kallisto_index   = getGenomeAttribute('kallisto')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                  } from './workflows/rnaseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_RNASEQ {

    main:

    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.ribo_database_manifest,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.gencode,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        params.skip_bbsplit,
        params.remove_ribo_rna ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { fai -> checkMaxContigSize(fai) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))

    // Bowtie2 rRNA index is built on-demand inside the fastq_remove_rrna subworkflow
    // rather than in PREPARE_GENOME, to avoid duplicating the rRNA FASTA preparation logic
    ch_bowtie2_index = channel.empty()

    RNASEQ (
        ch_samplesheet,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.rsem_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.kallisto_index,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.rrna_fastas,
        PREPARE_GENOME.out.sortmerna_index,
        ch_bowtie2_index,
        PREPARE_GENOME.out.splicesites
    )
    ch_versions = ch_versions.mix(RNASEQ.out.versions)

    emit:
    trim_status       = RNASEQ.out.trim_status       // channel: [id, boolean]
    map_status        = RNASEQ.out.map_status        // channel: [id, boolean]
    strand_status     = RNASEQ.out.strand_status     // channel: [id, boolean]
    multiqc_report    = RNASEQ.out.multiqc_report    // channel: /path/to/multiqc_report.html
    versions          = ch_versions                  // channel: [version1, version2, ...]
    collated_versions = RNASEQ.out.collated_versions // channel: /path/to/collated_versions.yml

    // Genome outputs
    fasta          = PREPARE_GENOME.out.fasta
    gtf            = PREPARE_GENOME.out.gtf
    fai            = PREPARE_GENOME.out.fai
    gene_bed       = PREPARE_GENOME.out.gene_bed
    transcript_fasta = PREPARE_GENOME.out.transcript_fasta
    chrom_sizes    = PREPARE_GENOME.out.chrom_sizes
    splicesites    = PREPARE_GENOME.out.splicesites
    star_index     = PREPARE_GENOME.out.star_index
    rsem_index     = PREPARE_GENOME.out.rsem_index
    hisat2_index   = PREPARE_GENOME.out.hisat2_index
    salmon_index   = PREPARE_GENOME.out.salmon_index
    kallisto_index = PREPARE_GENOME.out.kallisto_index
    bbsplit_index  = PREPARE_GENOME.out.bbsplit_index
    sortmerna_index = PREPARE_GENOME.out.sortmerna_index
    custom_fasta   = PREPARE_GENOME.out.custom_fasta
    custom_gtf     = PREPARE_GENOME.out.custom_gtf

    // QC and trimming outputs
    fastqc_raw_html    = RNASEQ.out.fastqc_raw_html
    fastqc_raw_zip     = RNASEQ.out.fastqc_raw_zip
    fastqc_trim_html   = RNASEQ.out.fastqc_trim_html
    fastqc_trim_zip    = RNASEQ.out.fastqc_trim_zip
    trim_html          = RNASEQ.out.trim_html
    trim_zip           = RNASEQ.out.trim_zip
    trim_log           = RNASEQ.out.trim_log
    trim_json          = RNASEQ.out.trim_json
    trim_unpaired      = RNASEQ.out.trim_unpaired
    umi_log            = RNASEQ.out.umi_log
    umi_reads          = RNASEQ.out.umi_reads
    umi_dedup          = RNASEQ.out.umi_dedup
    lint_log_raw       = RNASEQ.out.lint_log_raw
    lint_log_trimmed   = RNASEQ.out.lint_log_trimmed
    lint_log_bbsplit   = RNASEQ.out.lint_log_bbsplit
    lint_log_ribo      = RNASEQ.out.lint_log_ribo
    bbsplit_stats      = RNASEQ.out.bbsplit_stats
    sortmerna_log      = RNASEQ.out.sortmerna_log
    ribodetector_log   = RNASEQ.out.ribodetector_log
    seqkit_stats       = RNASEQ.out.seqkit_stats
    bowtie2_rrna_log   = RNASEQ.out.bowtie2_rrna_log
    bowtie2_rrna_index = RNASEQ.out.bowtie2_rrna_index
    seqkit_prefixed    = RNASEQ.out.seqkit_prefixed
    seqkit_converted   = RNASEQ.out.seqkit_converted

    // Alignment outputs
    star               = RNASEQ.out.star
    samtools           = RNASEQ.out.samtools
    transcriptome_bam  = RNASEQ.out.transcriptome_bam
    unaligned_sequences = RNASEQ.out.unaligned_sequences
    hisat2_summary     = RNASEQ.out.hisat2_summary
    samtools_bai       = RNASEQ.out.samtools_bai

    // MarkDuplicates outputs
    markdup            = RNASEQ.out.markdup

    // QC outputs
    preseq_txt         = RNASEQ.out.preseq_txt
    preseq_log         = RNASEQ.out.preseq_log
    qualimap_results   = RNASEQ.out.qualimap_results
    dupradar           = RNASEQ.out.dupradar

    // RSeQC outputs
    rseqc              = RNASEQ.out.rseqc

    // Contaminant screening outputs
    kraken_report      = RNASEQ.out.kraken_report
    bracken_txt        = RNASEQ.out.bracken_txt
    sylph_profile      = RNASEQ.out.sylph_profile
    sylphtax_output    = RNASEQ.out.sylphtax_output

    // Consolidated outputs
    stringtie_outputs     = RNASEQ.out.stringtie_outputs
    featurecounts_outputs = RNASEQ.out.featurecounts_outputs
    bigwig_outputs        = RNASEQ.out.bigwig_outputs
    pseudo             = RNASEQ.out.pseudo
    rsem               = RNASEQ.out.rsem
    star_salmon        = RNASEQ.out.star_salmon
    deseq2             = RNASEQ.out.deseq2
    pseudo_deseq2      = RNASEQ.out.pseudo_deseq2

    // MultiQC outputs
    multiqc_data           = RNASEQ.out.multiqc_data
    multiqc_plots          = RNASEQ.out.multiqc_plots
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        workflow.outputDir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        workflow.outputDir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RNASEQ.out.multiqc_report,
        NFCORE_RNASEQ.out.trim_status,
        NFCORE_RNASEQ.out.map_status,
        NFCORE_RNASEQ.out.strand_status
    )

    publish:
    // Genome outputs
    // Use ifEmpty([]) to handle empty collected channels (Nextflow bug #6745)
    fasta             = NFCORE_RNASEQ.out.fasta.ifEmpty([])
    gtf               = NFCORE_RNASEQ.out.gtf.ifEmpty([])
    fai               = NFCORE_RNASEQ.out.fai.ifEmpty([])
    gene_bed          = NFCORE_RNASEQ.out.gene_bed.ifEmpty([])
    transcript_fasta  = NFCORE_RNASEQ.out.transcript_fasta.ifEmpty([])
    chrom_sizes       = NFCORE_RNASEQ.out.chrom_sizes.ifEmpty([])
    splicesites       = NFCORE_RNASEQ.out.splicesites.ifEmpty([])
    star_index        = NFCORE_RNASEQ.out.star_index.ifEmpty([])
    rsem_index        = NFCORE_RNASEQ.out.rsem_index.ifEmpty([])
    hisat2_index      = NFCORE_RNASEQ.out.hisat2_index.ifEmpty([])
    salmon_index      = NFCORE_RNASEQ.out.salmon_index.ifEmpty([])
    kallisto_index    = NFCORE_RNASEQ.out.kallisto_index.ifEmpty([])
    bbsplit_index     = NFCORE_RNASEQ.out.bbsplit_index.ifEmpty([])
    sortmerna_index   = NFCORE_RNASEQ.out.sortmerna_index.ifEmpty([])
    custom_fasta      = NFCORE_RNASEQ.out.custom_fasta.ifEmpty([])
    custom_gtf        = NFCORE_RNASEQ.out.custom_gtf.ifEmpty([])

    // QC and trimming outputs
    fastqc_raw_html    = NFCORE_RNASEQ.out.fastqc_raw_html.ifEmpty([])
    fastqc_raw_zip     = NFCORE_RNASEQ.out.fastqc_raw_zip.ifEmpty([])
    fastqc_trim_html   = NFCORE_RNASEQ.out.fastqc_trim_html.ifEmpty([])
    fastqc_trim_zip    = NFCORE_RNASEQ.out.fastqc_trim_zip.ifEmpty([])
    trim_html          = NFCORE_RNASEQ.out.trim_html.ifEmpty([])
    trim_zip           = NFCORE_RNASEQ.out.trim_zip.ifEmpty([])
    trim_log           = NFCORE_RNASEQ.out.trim_log.ifEmpty([])
    trim_json          = NFCORE_RNASEQ.out.trim_json.ifEmpty([])
    trim_unpaired      = NFCORE_RNASEQ.out.trim_unpaired.ifEmpty([])
    umi_log            = NFCORE_RNASEQ.out.umi_log.ifEmpty([])
    umi_reads          = NFCORE_RNASEQ.out.umi_reads.ifEmpty([])
    umi_genomic_dedup_log        = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.genomic_dedup_log] }.ifEmpty([])
    umi_transcriptomic_dedup_log = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptomic_dedup_log] }.ifEmpty([])
    umi_prepare_for_rsem_log     = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.prepare_for_rsem_log] }.ifEmpty([])
    umi_transcriptome_dedup_bam      = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_dedup_bam] }.ifEmpty([])
    umi_transcriptome_sorted_bam     = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_sorted_bam] }.ifEmpty([])
    umi_transcriptome_sorted_bam_bai = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_sorted_bam_bai] }.ifEmpty([])
    umi_transcriptome_filtered_bam   = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_filtered_bam] }.ifEmpty([])
    umi_dedup_genome_stats    = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.genome_stats] }.ifEmpty([])
    umi_dedup_genome_flagstat = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.genome_flagstat] }.ifEmpty([])
    umi_dedup_genome_idxstats = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.genome_idxstats] }.ifEmpty([])
    umi_dedup_transcriptome_stats    = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_stats] }.ifEmpty([])
    umi_dedup_transcriptome_flagstat = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_flagstat] }.ifEmpty([])
    umi_dedup_transcriptome_idxstats = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.transcriptome_idxstats] }.ifEmpty([])
    umi_dedup_bam      = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.bam] }.ifEmpty([])
    umi_dedup_bai      = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.bai] }.ifEmpty([])
    umi_dedup_tsv_edit_distance    = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.tsv_edit_distance] }.ifEmpty([])
    umi_dedup_tsv_per_umi          = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.tsv_per_umi] }.ifEmpty([])
    umi_dedup_tsv_umi_per_position = NFCORE_RNASEQ.out.umi_dedup.map { r -> [r.meta, r.tsv_umi_per_position] }.ifEmpty([])
    lint_log_raw       = NFCORE_RNASEQ.out.lint_log_raw.ifEmpty([])
    lint_log_trimmed   = NFCORE_RNASEQ.out.lint_log_trimmed.ifEmpty([])
    lint_log_bbsplit   = NFCORE_RNASEQ.out.lint_log_bbsplit.ifEmpty([])
    lint_log_ribo      = NFCORE_RNASEQ.out.lint_log_ribo.ifEmpty([])
    bbsplit_stats      = NFCORE_RNASEQ.out.bbsplit_stats.ifEmpty([])
    sortmerna_log      = NFCORE_RNASEQ.out.sortmerna_log.ifEmpty([])
    ribodetector_log   = NFCORE_RNASEQ.out.ribodetector_log.ifEmpty([])
    seqkit_stats       = NFCORE_RNASEQ.out.seqkit_stats.ifEmpty([])
    bowtie2_rrna_log   = NFCORE_RNASEQ.out.bowtie2_rrna_log.ifEmpty([])
    bowtie2_rrna_index = NFCORE_RNASEQ.out.bowtie2_rrna_index.ifEmpty([])
    seqkit_prefixed    = NFCORE_RNASEQ.out.seqkit_prefixed.ifEmpty([])
    seqkit_converted   = NFCORE_RNASEQ.out.seqkit_converted.ifEmpty([])

    // Alignment outputs - extract from StarAlignResult and SamtoolsResult records
    star_log           = NFCORE_RNASEQ.out.star.map { r -> [r.meta, r.log_final] }.ifEmpty([])
    star_log_out       = NFCORE_RNASEQ.out.star.map { r -> [r.meta, r.log_out] }.ifEmpty([])
    star_log_progress  = NFCORE_RNASEQ.out.star.map { r -> [r.meta, r.log_progress] }.ifEmpty([])
    star_tab           = NFCORE_RNASEQ.out.star.map { r -> [r.meta, r.tab] }.ifEmpty([])
    star_bam           = NFCORE_RNASEQ.out.star.map { r -> [r.meta, r.bam] }.ifEmpty([])
    star_bai           = NFCORE_RNASEQ.out.samtools.map { r -> [r.meta, r.bai] }.ifEmpty([])
    sorted_bam_stats    = NFCORE_RNASEQ.out.samtools.map { r -> [r.meta, r.stats] }.ifEmpty([])
    sorted_bam_flagstat = NFCORE_RNASEQ.out.samtools.map { r -> [r.meta, r.flagstat] }.ifEmpty([])
    sorted_bam_idxstats = NFCORE_RNASEQ.out.samtools.map { r -> [r.meta, r.idxstats] }.ifEmpty([])
    transcriptome_bam  = NFCORE_RNASEQ.out.transcriptome_bam.ifEmpty([])
    unaligned_sequences = NFCORE_RNASEQ.out.unaligned_sequences.ifEmpty([])
    hisat2_summary     = NFCORE_RNASEQ.out.hisat2_summary.ifEmpty([])
    samtools_bai       = NFCORE_RNASEQ.out.samtools_bai.ifEmpty([])

    // MarkDuplicates outputs - extract from MarkDupResult record
    markdup_bam        = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.bam] }.ifEmpty([])
    markdup_bai        = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.bai] }.ifEmpty([])
    markdup_metrics    = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.metrics] }.ifEmpty([])
    markdup_stats      = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.stats] }.ifEmpty([])
    markdup_flagstat   = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.flagstat] }.ifEmpty([])
    markdup_idxstats   = NFCORE_RNASEQ.out.markdup.map { r -> [r.meta, r.idxstats] }.ifEmpty([])

    // QC outputs
    preseq_txt         = NFCORE_RNASEQ.out.preseq_txt.ifEmpty([])
    preseq_log         = NFCORE_RNASEQ.out.preseq_log.ifEmpty([])
    qualimap_results   = NFCORE_RNASEQ.out.qualimap_results.ifEmpty([])
    dupradar_scatter   = NFCORE_RNASEQ.out.dupradar.map { r -> [r.meta, r.scatter] }.ifEmpty([])
    dupradar_boxplot   = NFCORE_RNASEQ.out.dupradar.map { r -> [r.meta, r.boxplot] }.ifEmpty([])
    dupradar_histogram = NFCORE_RNASEQ.out.dupradar.map { r -> [r.meta, r.histogram] }.ifEmpty([])
    dupradar_gene_data = NFCORE_RNASEQ.out.dupradar.map { r -> [r.meta, r.gene_data] }.ifEmpty([])
    dupradar_intercept = NFCORE_RNASEQ.out.dupradar.map { r -> [r.meta, r.intercept] }.ifEmpty([])

    // RSeQC outputs - use .map to extract fields; null-safe ?. for nested records
    rseqc_bamstat              = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.bamstat] }.ifEmpty([])
    rseqc_inferexperiment      = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inferexperiment] }.ifEmpty([])
    rseqc_junctionannotation_bed         = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.bed] }.ifEmpty([])
    rseqc_junctionannotation_interact_bed = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.interact_bed] }.ifEmpty([])
    rseqc_junctionannotation_xls          = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.xls] }.ifEmpty([])
    rseqc_junctionannotation_log   = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.log] }.ifEmpty([])
    rseqc_junctionannotation_pdf        = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.pdf] }.ifEmpty([])
    rseqc_junctionannotation_events_pdf = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.events_pdf] }.ifEmpty([])
    rseqc_junctionannotation_r          = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junction_annotation?.rscript] }.ifEmpty([])
    rseqc_junctionsaturation_pdf   = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junctionsaturation_pdf] }.ifEmpty([])
    rseqc_junctionsaturation_r     = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.junctionsaturation_r] }.ifEmpty([])
    rseqc_readduplication_pos_xls  = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.read_duplication?.pos_xls] }.ifEmpty([])
    rseqc_readduplication_seq_xls  = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.read_duplication?.seq_xls] }.ifEmpty([])
    rseqc_readduplication_pdf      = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.read_duplication?.pdf] }.ifEmpty([])
    rseqc_readduplication_r        = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.read_duplication?.rscript] }.ifEmpty([])
    rseqc_readdistribution         = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.readdistribution] }.ifEmpty([])
    rseqc_innerdistance_txt        = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inner_distance?.freq] }.ifEmpty([])
    rseqc_innerdistance_distance   = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inner_distance?.distance] }.ifEmpty([])
    rseqc_innerdistance_mean       = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inner_distance?.mean] }.ifEmpty([])
    rseqc_innerdistance_pdf        = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inner_distance?.pdf] }.ifEmpty([])
    rseqc_innerdistance_r          = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.inner_distance?.rscript] }.ifEmpty([])
    rseqc_tin                      = NFCORE_RNASEQ.out.rseqc.map { r -> [r.meta, r.tin] }.ifEmpty([])

    // Contaminant screening outputs
    kraken_report      = NFCORE_RNASEQ.out.kraken_report.ifEmpty([])
    bracken_txt        = NFCORE_RNASEQ.out.bracken_txt.ifEmpty([])
    sylph_profile      = NFCORE_RNASEQ.out.sylph_profile.ifEmpty([])
    sylphtax_output    = NFCORE_RNASEQ.out.sylphtax_output.ifEmpty([])

    // Consolidated outputs - records publish all fields to their respective directories
    stringtie_outputs     = NFCORE_RNASEQ.out.stringtie_outputs.ifEmpty([])
    featurecounts_outputs = NFCORE_RNASEQ.out.featurecounts_outputs.ifEmpty([])
    bigwig_outputs        = NFCORE_RNASEQ.out.bigwig_outputs.ifEmpty([])
    pseudo                = NFCORE_RNASEQ.out.pseudo.ifEmpty([])
    rsem                  = NFCORE_RNASEQ.out.rsem.ifEmpty([])
    star_salmon           = NFCORE_RNASEQ.out.star_salmon.ifEmpty([])
    deseq2                = NFCORE_RNASEQ.out.deseq2.ifEmpty([])
    pseudo_deseq2         = NFCORE_RNASEQ.out.pseudo_deseq2.ifEmpty([])

    // MultiQC report
    multiqc_report         = NFCORE_RNASEQ.out.multiqc_report.ifEmpty([])
    multiqc_data           = NFCORE_RNASEQ.out.multiqc_data.ifEmpty([])
    multiqc_plots          = NFCORE_RNASEQ.out.multiqc_plots.ifEmpty([])

    // Pipeline info
    collated_versions      = NFCORE_RNASEQ.out.collated_versions.ifEmpty([])
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUT DEFINITIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    // Genome outputs - only publish if user wants references saved
    fasta            { enabled params.save_reference; path 'genome' }
    gtf              { enabled params.save_reference; path 'genome' }
    fai              { enabled params.save_reference; path 'genome' }
    gene_bed         { enabled params.save_reference; path 'genome' }
    transcript_fasta { enabled params.save_reference; path 'genome' }
    chrom_sizes      { enabled params.save_reference; path 'genome' }
    splicesites      { enabled params.save_reference; path 'genome' }
    star_index       { enabled params.save_reference; path 'genome/index/star' }
    rsem_index       { enabled params.save_reference; path 'genome/index/rsem' }
    hisat2_index     { enabled params.save_reference; path 'genome/index/hisat2' }
    salmon_index     { enabled params.save_reference; path 'genome/index/salmon' }
    kallisto_index   { enabled params.save_reference; path 'genome/index/kallisto' }
    bbsplit_index    { enabled params.save_reference; path 'genome/index/bbsplit' }
    sortmerna_index  { enabled params.save_reference; path 'genome/index/sortmerna' }
    custom_fasta     { path 'custom' }
    custom_gtf       { path 'custom' }

    // QC and trimming outputs
    fastqc_raw_html  { path 'fastqc/raw' }
    fastqc_raw_zip   { path 'fastqc/raw' }
    fastqc_trim_html { path 'fastqc/trim' }
    fastqc_trim_zip  { path 'fastqc/trim' }
    trim_html        { path { params.trimmer } }
    trim_zip         { path { params.trimmer } }
    trim_log         { path { params.trimmer == 'fastp' ? "${params.trimmer}/log" : params.trimmer } }
    trim_json        { path { params.trimmer } }
    trim_unpaired    { enabled params.save_trimmed; path { params.trimmer } }
    umi_log          { path 'umitools' }
    umi_reads        { enabled params.save_umi_intermeds; path 'umitools' }
    umi_genomic_dedup_log        { path { "${params.aligner}/${params.umi_dedup_tool}/genomic_dedup_log" } }
    umi_transcriptomic_dedup_log { path { "${params.aligner}/${params.umi_dedup_tool}/transcriptomic_dedup_log" } }
    umi_prepare_for_rsem_log     { path { "${params.aligner}/${params.umi_dedup_tool}/prepare_for_salmon_log" } }
    umi_transcriptome_dedup_bam      { path { params.aligner } }
    umi_transcriptome_sorted_bam     { path { params.aligner } }
    umi_transcriptome_sorted_bam_bai { path { params.aligner } }
    umi_transcriptome_filtered_bam   { path { params.aligner } }
    umi_dedup_genome_stats             { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_genome_flagstat          { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_genome_idxstats          { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_transcriptome_stats      { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_transcriptome_flagstat   { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_transcriptome_idxstats   { path { "${params.aligner}/samtools_stats" } }
    umi_dedup_bam    { path { params.aligner } }
    umi_dedup_bai    { path { params.aligner } }
    umi_dedup_tsv_edit_distance    { path { "${params.aligner}/${params.umi_dedup_tool}" } }
    umi_dedup_tsv_per_umi          { path { "${params.aligner}/${params.umi_dedup_tool}" } }
    umi_dedup_tsv_umi_per_position { path { "${params.aligner}/${params.umi_dedup_tool}" } }
    lint_log_raw     { path 'fq_lint/raw' }
    lint_log_trimmed { path 'fq_lint/trimmed' }
    lint_log_bbsplit { path 'fq_lint/bbsplit' }
    lint_log_ribo    { path { "fq_lint/${params.ribo_removal_tool ?: 'sortmerna'}" } }
    bbsplit_stats    { path 'bbsplit' }
    sortmerna_log    { path 'sortmerna' }
    ribodetector_log { path 'ribodetector' }
    seqkit_stats     { path 'ribodetector' }
    bowtie2_rrna_log { path 'bowtie2_rrna' }
    bowtie2_rrna_index { enabled params.save_reference; path 'genome/index/bowtie2_rrna' }
    seqkit_prefixed  { path 'seqkit' }
    seqkit_converted { path 'seqkit' }

    // Alignment outputs
    star_log            { path { "${params.aligner}/log" } }
    star_log_out        { path { "${params.aligner}/log" } }
    star_log_progress   { path { "${params.aligner}/log" } }
    star_tab            { path { "${params.aligner}/log" } }
    star_bam            { enabled params.save_align_intermeds; path { params.aligner } }
    star_bai            { enabled params.save_align_intermeds; path { params.aligner } }
    sorted_bam_stats    { path { "${params.aligner}/samtools_stats" } }
    sorted_bam_flagstat { path { "${params.aligner}/samtools_stats" } }
    sorted_bam_idxstats { path { "${params.aligner}/samtools_stats" } }
    transcriptome_bam   { enabled params.save_align_intermeds; path { params.aligner } }
    unaligned_sequences { enabled params.save_unaligned; path { "${params.aligner}/unmapped" } }
    hisat2_summary      { path { "${params.aligner}/log" } }
    samtools_bai        { path 'samtools' }  // Input BAM indices for BAM input mode

    // MarkDuplicates outputs
    markdup_bam      { path { params.aligner } }
    markdup_bai      { path { params.aligner } }
    markdup_metrics  { path { "${params.aligner}/picard_metrics" } }
    markdup_stats    { path { "${params.aligner}/samtools_stats" } }
    markdup_flagstat { path { "${params.aligner}/samtools_stats" } }
    markdup_idxstats { path { "${params.aligner}/samtools_stats" } }

    // QC outputs
    preseq_txt         { path { "${params.aligner}/preseq" } }
    preseq_log         { path { "${params.aligner}/preseq/log" } }
    qualimap_results   { path { "${params.aligner}/qualimap" } }
    dupradar_scatter   { path { "${params.aligner}/dupradar/scatter_plot" } }
    dupradar_boxplot   { path { "${params.aligner}/dupradar/box_plot" } }
    dupradar_histogram { path { "${params.aligner}/dupradar/histogram" } }
    dupradar_gene_data { path { "${params.aligner}/dupradar/gene_data" } }
    dupradar_intercept { path { "${params.aligner}/dupradar/intercepts_slope" } }

    // RSeQC outputs
    rseqc_bamstat                 { path { "${params.aligner}/rseqc/bam_stat" } }
    rseqc_inferexperiment         { path { "${params.aligner}/rseqc/infer_experiment" } }
    rseqc_junctionannotation_bed         { path { "${params.aligner}/rseqc/junction_annotation/bed" } }
    rseqc_junctionannotation_interact_bed { path { "${params.aligner}/rseqc/junction_annotation/bed" } }
    rseqc_junctionannotation_xls          { path { "${params.aligner}/rseqc/junction_annotation/xls" } }
    rseqc_junctionannotation_log  { path { "${params.aligner}/rseqc/junction_annotation/log" } }
    rseqc_junctionannotation_pdf        { path { "${params.aligner}/rseqc/junction_annotation/pdf" } }
    rseqc_junctionannotation_events_pdf { path { "${params.aligner}/rseqc/junction_annotation/pdf" } }
    rseqc_junctionannotation_r          { path { "${params.aligner}/rseqc/junction_annotation/rscript" } }
    rseqc_junctionsaturation_pdf  { path { "${params.aligner}/rseqc/junction_saturation/pdf" } }
    rseqc_junctionsaturation_r    { path { "${params.aligner}/rseqc/junction_saturation/rscript" } }
    rseqc_readduplication_pos_xls { path { "${params.aligner}/rseqc/read_duplication/xls" } }
    rseqc_readduplication_seq_xls { path { "${params.aligner}/rseqc/read_duplication/xls" } }
    rseqc_readduplication_pdf     { path { "${params.aligner}/rseqc/read_duplication/pdf" } }
    rseqc_readduplication_r       { path { "${params.aligner}/rseqc/read_duplication/rscript" } }
    rseqc_readdistribution        { path { "${params.aligner}/rseqc/read_distribution" } }
    rseqc_innerdistance_txt       { path { "${params.aligner}/rseqc/inner_distance/txt" } }
    rseqc_innerdistance_distance  { path { "${params.aligner}/rseqc/inner_distance/txt" } }
    rseqc_innerdistance_mean      { path { "${params.aligner}/rseqc/inner_distance/txt" } }
    rseqc_innerdistance_pdf       { path { "${params.aligner}/rseqc/inner_distance/pdf" } }
    rseqc_innerdistance_r         { path { "${params.aligner}/rseqc/inner_distance/rscript" } }
    rseqc_tin                     { path { "${params.aligner}/rseqc/tin" } }

    // Contaminant screening outputs
    kraken_report   { path { "${params.aligner}/contaminants/kraken2/kraken_reports" } }
    bracken_txt     { path { "${params.aligner}/contaminants/bracken" } }
    sylph_profile   { path { "${params.aligner}/contaminants/sylph" } }
    sylphtax_output { path { "${params.aligner}/contaminants/sylph" } }

    // Consolidated outputs - record channels publish all fields to one directory
    stringtie_outputs     { path { "${params.aligner}/stringtie" } }
    featurecounts_outputs { path { "${params.aligner}/featurecounts" } }
    bigwig_outputs        { path { "${params.aligner}/bigwig" } }
    pseudo                { path { params.pseudo_aligner } }
    star_salmon           { path 'star_salmon' }
    deseq2                { path { "${params.aligner}/deseq2_qc" } }
    pseudo_deseq2         { path { "${params.pseudo_aligner}/deseq2_qc" } }
    rsem                  { path 'star_rsem' }

    // MultiQC report
    multiqc_report { path { params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}" } }
    multiqc_data   { path { params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}" } }
    multiqc_plots  { path { params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}" } }

    // Pipeline info
    collated_versions { path 'pipeline_info' }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
