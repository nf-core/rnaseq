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
nextflow.preview.topic = true

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
    directory params.outdir
    mode params.publish_dir_mode

    'genome' {
        enabled params.save_reference
    }

    'genome/index' {
        enabled params.save_reference
    }

    'star_salmon/intermeds/' {
        path 'star_salmon'
        enabled params.save_align_intermeds || params.save_umi_intermeds
    }

    'star_salmon/samtools_stats/intermeds/' {
        path 'star_salmon/samtools_stats'
        enabled params.save_align_intermeds || params.save_umi_intermeds
    }

    'star_salmon/umitools/log/intermeds/' {
        path 'star_salmon/umitools/log'
        enabled params.save_align_intermeds || params.save_umi_intermeds
    }

    'bigwig' {
        path "${params.aligner}/bigwig"
    }

    // modules/local/dupradar
    'dupradar' {
        path "${params.aligner}/dupradar"
    }

    // modules/nf-core/bbmap/bbsplit
    'bbsplit/intermeds/' {
        path 'bbsplit'
        enabled params.save_bbsplit_reads
    }

    // modules/nf-core/cat/fastq
    'cat/fastq' {
        path 'fastq'
        enabled params.save_merged_fastq
    }

    // modules/nf-core/multiqc
    'multiqc' {
        path params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}"
    }

    // modules/nf-core/preseq/lcextrap
    'preseq' {
        path "${params.aligner}/preseq"
    }

    'preseq/log' {
        path "${params.aligner}/preseq/log"
    }

    // modules/nf-core/qualimap/rnaseq
    'qualimap' {
        path "${params.aligner}/qualimap"
    }

    // modules/nf-core/sortmerna
    'sortmerna/intermeds/' {
        path 'sortmerna'
        enabled params.save_non_ribo_reads
    }

    // modules/nf-core/stringtie/stringtie
    'stringtie' {
        path "${params.aligner}/stringtie"
    }

    // modules/nf-core/subread/featurecounts
    'featurecounts' {
        path "${params.aligner}/featurecounts"
    }

    // subworkflows/local/align_star
    // 'star_salmon/intermeds/' {
    //     path 'star_salmon'
    //     enabled params.save_align_intermeds
    // }

    'star_salmon/unmapped/' {
        enabled params.save_unaligned
    }

    // subworkflows/local/quantify_rsem
    'star_rsem/intermeds/' {
        path 'star_rsem'
        enabled params.save_align_intermeds
    }

    // subworkflows/nf-core/bam_markduplicates_picard
    'picard/' {
        path params.aligner
    }

    'picard/metrics/' {
        path "${params.aligner}/picard_metrics"
    }

    'picard/samtools_stats/' {
        path "${params.aligner}/samtools_stats"
    }

    // subworkflows/nf-core/bam_rseqc
    'rseqc/bam_stat/' {
        path "${params.aligner}/rseqc/bam_stat"
    }
    'rseqc/infer_experiment/' {
        path "${params.aligner}/rseqc/infer_experiment"
    }
    'rseqc/junction_annotation/' {
        path "${params.aligner}/rseqc/junction_annotation"
    }
    'rseqc/junction_saturation/' {
        path "${params.aligner}/rseqc/junction_saturation"
    }
    'rseqc/read_duplication/' {
        path "${params.aligner}/rseqc/read_duplication"
    }
    'rseqc/read_distribution/' {
        path "${params.aligner}/rseqc/read_distribution"
    }
    'rseqc/inner_distance/' {
        path "${params.aligner}/rseqc/inner_distance"
    }
    'rseqc/tin/' {
        path "${params.aligner}/rseqc/tin"
    }

    // subworkflows/nf-core/fastq_align_hisat2
    'hisat2/intermeds/' {
        path 'hisat2'
        enabled params.save_align_intermeds
    }

    'hisat2/unmapped/' {
        enabled params.save_unaligned
    }

    // subworkflows/nf-core/fastq_fastqc_umitools_fastp
    // subworkflows/nf-core/fastq_fastqc_umitools_trimgalore
    'fastp/intermeds/' {
        path 'fastp'
        enabled params.save_trimmed
    }

    'trimgalore/intermeds/' {
        path 'trimgalore'
        enabled params.save_trimmed
    }

    'umitools/intermeds/' {
        path 'umitools'
        enabled params.save_umi_intermeds
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
