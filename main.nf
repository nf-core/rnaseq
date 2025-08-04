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

nextflow.preview.output = true

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

include { RNASEQ                  } from './workflows/rnaseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

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
        params.outdir
    )

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
        !params.remove_ribo_rna,
        params.skip_alignment,
        params.skip_pseudo_alignment
    )

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { checkMaxContigSize(it) }
    }

    ch_genome = Channel.empty().mix(
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gff,
        PREPARE_GENOME.out.add_fasta,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
    )

    ch_genome_index = Channel.empty().mix(
        PREPARE_GENOME.out.splicesites,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.rsem_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.kallisto_index,
    )

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))
    RNASEQ (
        ch_samplesheet,
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
        PREPARE_GENOME.out.splicesites
    )

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
        RNASEQ.out.multiqc_report,
        RNASEQ.out.trim_status,
        RNASEQ.out.map_status,
        RNASEQ.out.strand_status
    )

    publish:
    genome = ch_genome
    genome_index = ch_genome_index
    star_salmon = RNASEQ.out.star_salmon
    star_salmon_deseq_qc = RNASEQ.out.star_salmon_deseq_qc
    star_rsem = RNASEQ.out.star_rsem
    star_rsem_deseq_qc = RNASEQ.out.star_rsem_deseq_qc
    hisat2 = RNASEQ.out.hisat2
    multiqc_report = RNASEQ.out.multiqc_report
    multiqc_data = RNASEQ.out.multiqc_data
    multiqc_plots = RNASEQ.out.multiqc_plots
}

output {
    genome {
        enabled params.save_reference
        path 'genome'
    }

    genome_index {
        enabled params.save_reference
        path 'genome/index'
    }

    star_salmon {
        path 'star_salmon'
    }

    star_salmon_deseq_qc {
        path 'star_salmon/deseq2_qc'
    }

    star_rsem {
        path 'star_rsem'
    }

    star_rsem_deseq_qc {
        path 'star_rsem/deseq2_qc'
    }

    hisat2 {
        path 'hisat2'
    }

    multiqc_report {
        path params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}"
    }

    multiqc_data {
        path params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}"
    }

    multiqc_plots {
        path params.skip_alignment ? 'multiqc' : "multiqc/${params.aligner}"
    }
}

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
