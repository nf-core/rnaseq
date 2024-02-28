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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                  } from './workflows/rnaseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.gtf              = getGenomeAttribute('gtf')
params.gff              = getGenomeAttribute('gff')
params.gene_bed         = getGenomeAttribute('bed12')
params.bbsplit_index    = getGenomeAttribute('bbsplit')
params.star_index       = getGenomeAttribute('star')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.rsem_index       = getGenomeAttribute('rsem')
params.salmon_index     = getGenomeAttribute('salmon')
params.kallisto_index   = getGenomeAttribute('kallisto')

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

    ch_versions = Channel.empty()

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
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.gencode,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        params.skip_bbsplit,
        params.skip_alignment,
        params.skip_pseudo_alignment
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { checkMaxContigSize(it) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))
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
        PREPARE_GENOME.out.splicesites
    )
    ch_versions = ch_versions.mix(RNASEQ.out.versions)

    emit:
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions               // channel: [version1, version2, ...]
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
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
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
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RNASEQ.out.multiqc_report
    )
}

output {
    path(params.outdir) {
        path('genome', enabled: params.save_reference) {
            select 'GUNZIP_.*'
            select 'MAKE_TRANSCRIPTS_FASTA'
            select 'GFFREAD'
            select 'GTF2BED'
            select 'CAT_ADDITIONAL_FASTA'
            select 'PREPROCESS_TRANSCRIPTS_FASTA_GENCODE'
            select 'GTF_FILTER'
            select 'CUSTOM_GETCHROMSIZES'
        }

        path('genome/index', enabled: params.save_reference) {
            select 'UNTAR_.*'
            select 'STAR_GENOMEGENERATE'
            select 'STAR_GENOMEGENERATE_IGENOMES'
            select 'HISAT2_BUILD'
            select 'HISAT2_EXTRACTSPLICESITES'
            select 'SALMON_INDEX'
            select 'KALLISTO_INDEX'
            select 'RSEM_PREPAREREFERENCE_GENOME'
            select 'PREPARE_GENOME:BBMAP_BBSPLIT'
        }

        path(params.aligner) {
            select '.*:QUANTIFY_STAR_SALMON:SALMON_QUANT'
            select '.*:QUANTIFY_STAR_SALMON:CUSTOM_TX2GENE'
            select '.*:QUANTIFY_STAR_SALMON:TXIMETA_TXIMPORT'
            select '.*:QUANTIFY_STAR_SALMON:SE_.*'

            path('umitools') {
                select '.*:BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME:UMITOOLS_DEDUP', pattern: '*.tsv'
            }

            path('samtools_stats') {
                select '.*:BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME:BAM_STATS_SAMTOOLS:.*', pattern: '*.{stats,flagstat,idxstats}'
                select 'NFCORE_RNASEQ:RNASEQ:.*:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:.*', pattern: '*.{stats,flagstat,idxstats}'
            }
        }

        path(params.aligner, enabled: params.save_align_intermeds || params.save_umi_intermeds) {
            select 'NFCORE_RNASEQ:RNASEQ:SAMTOOLS_SORT', pattern: '*.bam'
            select 'NFCORE_RNASEQ:RNASEQ:UMITOOLS_PREPAREFORSALMON', pattern: '*.bam'
            path('umitools/log') {
                select 'NFCORE_RNASEQ:RNASEQ:UMITOOLS_PREPAREFORSALMON', pattern: '*.log'
            }

            select 'NFCORE_RNASEQ:RNASEQ:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT', pattern: '*.bam'
            select 'NFCORE_RNASEQ:RNASEQ:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX', pattern: '*.bai'
            path('samtools_stats') {
                select 'NFCORE_RNASEQ:RNASEQ:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:.*', pattern: '*.{stats,flagstat,idxstats}'
            }

            select '.*:BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME:UMITOOLS_DEDUP', pattern: '*.bam'
            select '.*:BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME:SAMTOOLS_INDEX', pattern: '*.bai'
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
