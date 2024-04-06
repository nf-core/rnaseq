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
    directory params.outdir, mode: params.publish_dir_mode

    //
    // Genome preparation
    //
    'genome' {
        enabled params.save_reference
        from 'genome'
    }

    'genome/index' {
        enabled params.save_reference
        from 'genome-index'
    }

    //
    // Alignment
    //
    "${params.aligner}" {
        from 'align'

        'samtools_stats' {
            from 'align-samtools-stats'
        }

        'umitools' {
            from 'align-umitools'
        }
    }

    "${params.aligner}" {
        enabled params.save_align_intermeds || params.save_umi_intermeds
        from 'align-intermeds'

        'samtools_stats' {
            from 'align-intermeds-samtools-stats'
        }
        'umitools/log' {
            from 'align-intermeds-umitools-log'
        }
    }

    //
    // bigWig coverage
    //
    "${params.aligner}/bigwig" {
        from 'align-bigwig'
    }

    //
    // DESeq2 QC
    //
    "${params.aligner}/deseq2_qc" {
        from 'align-deseq2'
    }

    //
    // Pseudo-alignment
    //
    "${params.pseudo_aligner}" {
        from 'pseudo-align'
        // TODO: process 'QUANTIFY_PSEUDO_ALIGNMENT:CUSTOM_TX2GENE'

        'deseq2_qc' {
            from 'pseudo-align-deseq2'
        }
    }

    // modules/local/dupradar
    "${params.aligner}/dupradar" {
        'scatter_plot' {
            from 'align-dupradar', pattern: "*Dens.pdf"
        }
        'box_plot' {
            from 'align-dupradar', pattern: "*Boxplot.pdf"
        }
        'histogram' {
            from 'align-dupradar', pattern: "*Hist.pdf"
        }
        'gene_data' {
            from 'align-dupradar', pattern: "*Matrix.txt"
        }
        'intercepts_slope' {
            from 'align-dupradar', pattern: "*slope.txt"
        }
    }

    // modules/nf-core/bbmap/bbsplit
    'bbsplit' {
        from 'bbsplit'
        from 'bbsplit-intermeds', enabled: params.save_bbsplit_reads
    }

    // modules/nf-core/cat/fastq
    'fastq' {
        enabled params.save_merged_fastq
        from 'merged-fastq'
    }

    // modules/nf-core/multiqc
    def multiqc_path = params.skip_alignment
        ? 'multiqc'
        : "multiqc/${params.aligner}"

    "${multiqc_path}" {
        from 'multiqc'
    }

    // modules/nf-core/preseq/lcextrap
    "${params.aligner}" {
        'preseq' {
            from 'preseq'
        }
        'preseq/log' {
            from 'preseq-log'
        }
    }

    // modules/nf-core/qualimap/rnaseq
    "${params.aligner}/qualimap" {
        from 'qualimap'
    }

    // modules/nf-core/sortmerna
    'sortmerna' {
        from 'sortmerna'
        from 'sortmerna-intermeds', enabled: params.save_non_ribo_reads
    }

    // modules/nf-core/stringtie/stringtie
    "${params.aligner}/stringtie" {
        from 'align-stringtie'
    }

    // modules/nf-core/subread/featurecounts
    "${params.aligner}/featurecounts" {
        from 'align-featurecounts'
    }

    // subworkflows/local/align_star
    "${params.aligner}" {
        'log' {
            from 'align-star-log'
        }

        from 'align-star-intermeds'

        'unmapped' {
            from 'align-star-unaligned'
        }
    }

    // subworkflows/local/quantify_rsem
    "${params.aligner}" {
        from 'quantify-rsem'
        from 'quantify-rsem-intermeds', enabled: params.save_align_intermeds

        'log' {
            from 'quantify-rsem-log'
        }
    }

    // subworkflows/nf-core/bam_markduplicates_picard
    "${params.aligner}" {
        from 'align-picard'
        'picard_metrics' {
            from 'align-picard-metrics'
        }
        'samtools_stats' {
            from 'align-picard-stats'
        }
    }

    // subworkflows/nf-core/bam_rseqc
    "${params.aligner}/rseqc" {
        'bam_stat' {
            from 'align-rseqc-bamstat'
        }
        'infer_experiment' {
            from 'align-rseqc-inferexperiment'
        }
        'junction_annotation' {
            from 'align-rseqc-junctionannotation'
        }
        'junction_saturation' {
            from 'align-rseqc-junctionsaturation'
        }
        'read_duplication' {
            from 'align-rseqc-readduplication'
        }
        'read_distribution' {
            from 'align-rseqc-readdistribution'
        }
        'inner_distance' {
            from 'align-rseqc-innerdistance'
        }
        'tin' {
            from 'align-rseqc-tin'
        }
    }

    // subworkflows/nf-core/fastq_align_hisat2
    "${params.aligner}" {
        'log' {
            from 'align-hisat2-log'
        }

        from 'align-hisat2-intermeds', enabled: params.save_align_intermeds

        'unmapped' {
            from 'align-hisat2-unaligned', enabled: params.save_unaligned
        }
    }

    // subworkflows/nf-core/fastq_fastqc_umitools_fastp
    // subworkflows/nf-core/fastq_fastqc_umitools_trimgalore
    "${params.trimmer}" {
        'fastqc' {
            from 'trim-fastqc'
        }
        from 'trim'
        'log' {
            from 'trim-log'
        }
        from 'trim-intermeds', enabled: params.save_trimmed
    }

    'umitools' {
        from 'umitools'
        from 'umitools-intermeds', enabled: params.save_umi_intermeds
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
