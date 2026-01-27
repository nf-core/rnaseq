/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../../modules/local/deseq2_qc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../../modules/local/multiqc_custom_biotype'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ALIGN_STAR                            } from '../../subworkflows/local/align_star'
include { ALIGN_HISAT2                          } from '../../subworkflows/local/align_hisat2'
include { QUANTIFY_RSEM                         } from '../../subworkflows/local/quantify_rsem'

include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList             } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { biotypeInGtf                   } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getInferexperimentStrandedness } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { methodsDescriptionText         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { mapBamToPublishedPath          } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { DUPRADAR                   } from '../../modules/nf-core/dupradar'
include { PRESEQ_LCEXTRAP            } from '../../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ            } from '../../modules/nf-core/qualimap/rnaseq'
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie'
include { SUBREAD_FEATURECOUNTS      } from '../../modules/nf-core/subread/featurecounts'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { SYLPH_PROFILE              } from '../../modules/nf-core/sylph/profile/main'
include { SYLPHTAX_TAXPROF           } from '../../modules/nf-core/sylphtax/taxprof/main'
include { MULTIQC                    } from '../../modules/nf-core/multiqc'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QUALIMAP              } from '../../modules/nf-core/samtools/sort'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../../subworkflows/nf-core/bam_rseqc'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_STAR_SALMON } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet       // channel: path(sample_sheet.csv)
    ch_versions          // channel: [ path(versions.yml) ]
    ch_fasta             // channel: path(genome.fasta)
    ch_gtf               // channel: path(genome.gtf)
    ch_fai               // channel: path(genome.fai)
    ch_chrom_sizes       // channel: path(genome.sizes)
    ch_gene_bed          // channel: path(gene.bed)
    ch_transcript_fasta  // channel: path(transcript.fasta)
    ch_star_index        // channel: path(star/index/)
    ch_rsem_index        // channel: path(rsem/index/)
    ch_hisat2_index      // channel: path(hisat2/index/)
    ch_salmon_index      // channel: path(salmon/index/)
    ch_kallisto_index    // channel: [ meta, path(kallisto/index/) ]
    ch_bbsplit_index     // channel: path(bbsplit/index/)
    ch_ribo_db           // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index   // channel: path(sortmerna/index/)
    ch_bowtie2_index     // channel: path(bowtie2/index/) for rRNA removal
    ch_splicesites       // channel: path(genome.splicesites.txt)

    main:

    // Header files for MultiQC
    def ch_pca_header_multiqc        = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
    def sample_status_header_multiqc = file("$projectDir/workflows/rnaseq/assets/multiqc/sample_status_header.txt", checkIfExists: true)
    def ch_clustering_header_multiqc = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
    def ch_biotypes_header_multiqc   = file("$projectDir/workflows/rnaseq/assets/multiqc/biotypes_header.txt", checkIfExists: true)
    def ch_dummy_file                = ch_pca_header_multiqc

    ch_multiqc_files = channel.empty()
    ch_trim_status = channel.empty()
    ch_map_status = channel.empty()
    ch_strand_status = channel.empty()
    ch_percent_mapped = channel.empty()

    // Channel initializations for workflow outputs
    ch_stringtie_gtf          = channel.empty()
    ch_stringtie_coverage     = channel.empty()
    ch_stringtie_abundance    = channel.empty()
    ch_stringtie_ballgown     = channel.empty()
    ch_featurecounts_counts   = channel.empty()
    ch_featurecounts_summary  = channel.empty()
    ch_biotype_counts         = channel.empty()
    ch_bigwig_forward         = channel.empty()
    ch_bigwig_reverse         = channel.empty()
    ch_preseq_txt             = channel.empty()
    ch_preseq_log             = channel.empty()
    ch_markdup_bam            = channel.empty()
    ch_markdup_bai            = channel.empty()
    ch_markdup_metrics        = channel.empty()
    ch_markdup_stats          = channel.empty()
    ch_markdup_flagstat       = channel.empty()
    ch_markdup_idxstats       = channel.empty()
    ch_qualimap_results       = channel.empty()
    ch_dupradar_scatter       = channel.empty()
    ch_dupradar_boxplot       = channel.empty()
    ch_dupradar_histogram     = channel.empty()
    ch_dupradar_gene_data     = channel.empty()
    ch_dupradar_intercept     = channel.empty()
    ch_rseqc_bamstat          = channel.empty()
    ch_rseqc_inferexperiment  = channel.empty()
    ch_rseqc_junctionannotation_bed   = channel.empty()
    ch_rseqc_junctionannotation_xls        = channel.empty()
    ch_rseqc_junctionannotation_log        = channel.empty()
    ch_rseqc_junctionannotation_pdf        = channel.empty()
    ch_rseqc_junctionannotation_events_pdf = channel.empty()
    ch_rseqc_junctionannotation_r          = channel.empty()
    ch_rseqc_junctionsaturation_pdf   = channel.empty()
    ch_rseqc_junctionsaturation_r     = channel.empty()
    ch_rseqc_readduplication_pos_xls  = channel.empty()
    ch_rseqc_readduplication_seq_xls  = channel.empty()
    ch_rseqc_readduplication_pdf      = channel.empty()
    ch_rseqc_readduplication_r        = channel.empty()
    ch_rseqc_readdistribution                = channel.empty()
    ch_rseqc_innerdistance_txt               = channel.empty()
    ch_rseqc_innerdistance_distance          = channel.empty()
    ch_rseqc_innerdistance_mean              = channel.empty()
    ch_rseqc_innerdistance_pdf               = channel.empty()
    ch_rseqc_innerdistance_r                 = channel.empty()
    ch_rseqc_tin                             = channel.empty()
    ch_rseqc_junctionannotation_interact_bed = channel.empty()
    ch_kraken_report          = channel.empty()
    ch_bracken_txt            = channel.empty()
    ch_sylph_profile          = channel.empty()
    ch_sylphtax_output        = channel.empty()
    ch_pseudo_quant           = channel.empty()
    ch_pseudo_tx2gene         = channel.empty()
    ch_pseudo_counts_gene     = channel.empty()
    ch_pseudo_counts_gene_length_scaled = channel.empty()
    ch_pseudo_counts_gene_scaled = channel.empty()
    ch_pseudo_counts_transcript = channel.empty()
    ch_pseudo_lengths_gene    = channel.empty()
    ch_pseudo_lengths_transcript = channel.empty()
    ch_pseudo_tpm_gene        = channel.empty()
    ch_pseudo_tpm_transcript  = channel.empty()
    ch_pseudo_merged_gene_rds = channel.empty()
    ch_pseudo_merged_transcript_rds = channel.empty()
    ch_rsem_stat              = channel.empty()
    ch_rsem_logs              = channel.empty()
    ch_rsem_counts_gene       = channel.empty()
    ch_rsem_counts_transcript = channel.empty()
    ch_rsem_tpm_gene          = channel.empty()
    ch_rsem_tpm_transcript    = channel.empty()
    ch_rsem_merged_counts_gene       = channel.empty()
    ch_rsem_merged_counts_transcript = channel.empty()
    ch_rsem_merged_genes_long        = channel.empty()
    ch_rsem_merged_isoforms_long     = channel.empty()
    ch_star_salmon_quant      = channel.empty()
    ch_star_salmon_tx2gene    = channel.empty()
    ch_star_salmon_counts_gene = channel.empty()
    ch_star_salmon_counts_gene_length_scaled = channel.empty()
    ch_star_salmon_counts_gene_scaled = channel.empty()
    ch_star_salmon_counts_transcript = channel.empty()
    ch_star_salmon_lengths_gene = channel.empty()
    ch_star_salmon_lengths_transcript = channel.empty()
    ch_star_salmon_tpm_gene   = channel.empty()
    ch_star_salmon_tpm_transcript = channel.empty()
    ch_star_salmon_merged_gene_rds       = channel.empty()
    ch_star_salmon_merged_transcript_rds = channel.empty()
    ch_deseq2_pca             = channel.empty()
    ch_deseq2_dists           = channel.empty()
    ch_deseq2_pdf             = channel.empty()
    ch_deseq2_rdata           = channel.empty()
    ch_deseq2_pca_txt         = channel.empty()
    ch_deseq2_dists_txt       = channel.empty()
    ch_deseq2_log             = channel.empty()
    ch_deseq2_size_factors    = channel.empty()
    // Separate channels for pseudo-aligner DESeq2 QC (goes to different output path)
    ch_pseudo_deseq2_pca             = channel.empty()
    ch_pseudo_deseq2_dists           = channel.empty()
    ch_pseudo_deseq2_pdf             = channel.empty()
    ch_pseudo_deseq2_rdata           = channel.empty()
    ch_pseudo_deseq2_pca_txt         = channel.empty()
    ch_pseudo_deseq2_dists_txt       = channel.empty()
    ch_pseudo_deseq2_log             = channel.empty()
    ch_pseudo_deseq2_size_factors    = channel.empty()
    ch_hisat2_summary         = channel.empty()
    ch_star_bam               = channel.empty()
    ch_star_bai               = channel.empty()
    ch_sorted_bam_stats       = channel.empty()
    ch_sorted_bam_flagstat    = channel.empty()
    ch_sorted_bam_idxstats    = channel.empty()
    ch_transcriptome_bam_out  = channel.empty()
    ch_umi_genomic_dedup_log        = channel.empty()
    ch_umi_transcriptomic_dedup_log = channel.empty()
    ch_umi_prepare_for_rsem_log     = channel.empty()
    ch_umi_transcriptome_dedup_bam      = channel.empty()
    ch_umi_transcriptome_sorted_bam     = channel.empty()
    ch_umi_transcriptome_sorted_bam_bai = channel.empty()
    ch_umi_transcriptome_filtered_bam   = channel.empty()
    ch_umi_dedup_stats        = channel.empty()
    ch_umi_dedup_bam          = channel.empty()
    ch_umi_dedup_bai          = channel.empty()
    ch_umi_dedup_flagstat     = channel.empty()
    ch_umi_dedup_idxstats     = channel.empty()
    ch_umi_dedup_tsv_edit_distance    = channel.empty()
    ch_umi_dedup_tsv_per_umi          = channel.empty()
    ch_umi_dedup_tsv_umi_per_position = channel.empty()
    ch_samtools_bai           = channel.empty()

    //
    // Collect versions from topic channel (for modules that emit versions via topics)
    //
    def topic_versions = channel.topic('versions')
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    //
    // Create channel from input file provided through params.input
    //
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2, genome_bam, transcriptome_bam ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ], genome_bam, transcriptome_bam ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ], genome_bam, transcriptome_bam ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            checkSamplesAfterGrouping(samplesheet)
        }
        .branch {
            meta, reads, genome_bam, transcriptome_bam ->
                bam: params.skip_alignment && (genome_bam || transcriptome_bam)
                    return [ meta, genome_bam, transcriptome_bam ]
                fastq: reads.size() > 0 && reads[0]
                    return [ meta.findAll { key, _value -> key != 'percent_mapped' }, reads ]
        }
        .set { ch_input_branched }

    // Get inputs for FASTQ and BAM processing paths

    ch_fastq = ch_input_branched.fastq
    ch_genome_bam = ch_input_branched.bam.map { meta, genome_bam, _transcriptome_bam -> [ meta, genome_bam ] }.distinct()
    ch_transcriptome_bam = ch_input_branched.bam.map { meta, _genome_bam, transcriptome_bam -> [ meta, transcriptome_bam ] }.distinct()

    // Derive mapping percentages if supplied with input

    ch_percent_mapped = ch_input_branched.bam
        .filter{ meta, _genome_bam, _transcriptome_bam -> meta.percent_mapped }
        .map { meta, _genome_bam, _transcriptome_bam -> [ meta, meta.percent_mapped ] }

    // Index pre-aligned input BAM files
    SAMTOOLS_INDEX (
        ch_genome_bam
    )
    ch_genome_bam_index = params.bam_csi_index ? SAMTOOLS_INDEX.out.csi : SAMTOOLS_INDEX.out.bai
    ch_samtools_bai     = params.bam_csi_index ? SAMTOOLS_INDEX.out.csi : SAMTOOLS_INDEX.out.bai // For publishing input BAM indices
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    // Determine if we need to build rRNA removal indexes
    def make_sortmerna_index = !params.sortmerna_index && params.remove_ribo_rna && params.ribo_removal_tool == 'sortmerna'
    def make_bowtie2_index   = params.remove_ribo_rna && params.ribo_removal_tool == 'bowtie2'

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,                                   // ch_reads
        ch_fasta,                                   // ch_fasta
        ch_transcript_fasta,                        // ch_transcript_fasta
        ch_gtf,                                     // ch_gtf
        ch_salmon_index,                            // ch_salmon_index
        ch_sortmerna_index,                         // ch_sortmerna_index
        ch_bowtie2_index,                           // ch_bowtie2_index
        ch_bbsplit_index,                           // ch_bbsplit_index
        ch_ribo_db,                                 // ch_rrna_fastas
        params.skip_bbsplit || !params.fasta,       // skip_bbsplit
        params.skip_fastqc || params.skip_qc,       // skip_fastqc
        params.skip_trimming,                       // skip_trimming
        params.skip_umi_extract,                    // skip_umi_extract
        params.skip_linting,                        // skip_linting
        !salmon_index_available,                    // make_salmon_index
        make_sortmerna_index,                       // make_sortmerna_index
        make_bowtie2_index,                         // make_bowtie2_index
        params.trimmer,                             // trimmer
        params.min_trimmed_reads,                   // min_trimmed_reads
        params.save_trimmed,                        // save_trimmed
        false,                                      // fastp_merge
        params.remove_ribo_rna,                     // remove_ribo_rna
        params.ribo_removal_tool,                   // ribo_removal_tool
        params.with_umi,                            // with_umi
        params.umi_discard_read,                    // umi_discard_read
        params.stranded_threshold,                  // stranded_threshold
        params.unstranded_threshold                 // unstranded_threshold
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    // Capture individual outputs for workflow outputs
    ch_fastqc_raw_html    = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.fastqc_raw_html
    ch_fastqc_raw_zip     = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.fastqc_raw_zip
    ch_fastqc_trim_html   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.fastqc_trim_html
    ch_fastqc_trim_zip    = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.fastqc_trim_zip
    ch_trim_html          = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_html
    ch_trim_zip           = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_zip
    ch_trim_log           = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_log
    ch_trim_json          = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_json
    ch_trim_unpaired      = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_unpaired
    ch_umi_log            = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.umi_log
    ch_umi_reads          = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.umi_reads
    ch_bbsplit_stats      = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.bbsplit_stats
    ch_sortmerna_log      = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.sortmerna_log
    ch_ribodetector_log   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.ribodetector_log
    ch_seqkit_stats       = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.seqkit_stats
    ch_bowtie2_rrna_log   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.bowtie2_log
    ch_bowtie2_rrna_index = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.bowtie2_index
    ch_seqkit_prefixed    = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.seqkit_prefixed
    ch_seqkit_converted   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.seqkit_converted
    ch_lint_log_raw       = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.lint_log_raw
    ch_lint_log_trimmed   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.lint_log_trimmed
    ch_lint_log_bbsplit   = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.lint_log_bbsplit
    ch_lint_log_ribo      = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.lint_log_ribo

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_star_log            = channel.empty()
    ch_star_log_out        = channel.empty()
    ch_star_log_progress   = channel.empty()
    ch_star_tab            = channel.empty()
    ch_unaligned_sequences = channel.empty()

    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem')) {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        ALIGN_STAR (
            ch_strand_inferred_filtered_fastq,
            ch_star_index.map { item -> [ [:], item ] },
            ch_gtf.map { item -> [ [:], item ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_fasta.map { item -> [ [:], item ] },
            params.use_sentieon_star,
            params.with_umi,
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            params.bam_csi_index,
            params.skip_markduplicates,
            ch_transcript_fasta.map { item -> [ [:], item ] },
            ch_genome_bam,
            ch_genome_bam_index,
            ch_transcriptome_bam
        )

        ch_genome_bam            = ALIGN_STAR.out.bam
        ch_genome_bam_index      = ALIGN_STAR.out.bai
        ch_transcriptome_bam     = ALIGN_STAR.out.bam_transcript
        ch_transcriptome_bam_out = ALIGN_STAR.out.bam_transcript
        ch_percent_mapped        = ch_percent_mapped.mix(ALIGN_STAR.out.percent_mapped)
        ch_unprocessed_bams              = ALIGN_STAR.out.orig_bam.join(ALIGN_STAR.out.bam_transcript)
        ch_star_log                      = ALIGN_STAR.out.log_final
        ch_star_log_out                  = ALIGN_STAR.out.log_out
        ch_star_log_progress             = ALIGN_STAR.out.log_progress
        ch_star_tab                      = ALIGN_STAR.out.tab
        ch_star_bam                      = ALIGN_STAR.out.bam
        ch_star_bai                      = ALIGN_STAR.out.bai
        ch_sorted_bam_stats              = ALIGN_STAR.out.stats
        ch_sorted_bam_flagstat           = ALIGN_STAR.out.flagstat
        ch_sorted_bam_idxstats           = ALIGN_STAR.out.idxstats
        ch_unaligned_sequences           = ALIGN_STAR.out.fastq
        ch_umi_genomic_dedup_log        = ch_umi_genomic_dedup_log.mix(ALIGN_STAR.out.umi_genomic_dedup_log)
        ch_umi_transcriptomic_dedup_log = ch_umi_transcriptomic_dedup_log.mix(ALIGN_STAR.out.umi_transcriptomic_dedup_log)
        ch_umi_prepare_for_rsem_log     = ch_umi_prepare_for_rsem_log.mix(ALIGN_STAR.out.umi_prepare_for_rsem_log)
        ch_umi_transcriptome_dedup_bam      = ch_umi_transcriptome_dedup_bam.mix(ALIGN_STAR.out.umi_transcriptome_dedup_bam)
        ch_umi_transcriptome_sorted_bam     = ch_umi_transcriptome_sorted_bam.mix(ALIGN_STAR.out.umi_transcriptome_sorted_bam)
        ch_umi_transcriptome_sorted_bam_bai = ch_umi_transcriptome_sorted_bam_bai.mix(ALIGN_STAR.out.umi_transcriptome_sorted_bam_bai)
        ch_umi_transcriptome_filtered_bam   = ch_umi_transcriptome_filtered_bam.mix(ALIGN_STAR.out.umi_transcriptome_filtered_bam)
        ch_umi_dedup_stats    = ch_umi_dedup_stats.mix(ALIGN_STAR.out.umi_dedup_stats)
        ch_umi_dedup_bam      = ch_umi_dedup_bam.mix(ALIGN_STAR.out.umi_dedup_bam)
        ch_umi_dedup_bai      = ch_umi_dedup_bai.mix(ALIGN_STAR.out.umi_dedup_bai)
        ch_umi_dedup_flagstat = ch_umi_dedup_flagstat.mix(ALIGN_STAR.out.umi_dedup_flagstat)
        ch_umi_dedup_idxstats = ch_umi_dedup_idxstats.mix(ALIGN_STAR.out.umi_dedup_idxstats)
        ch_umi_dedup_tsv_edit_distance    = ch_umi_dedup_tsv_edit_distance.mix(ALIGN_STAR.out.umi_dedup_tsv_edit_distance)
        ch_umi_dedup_tsv_per_umi          = ch_umi_dedup_tsv_per_umi.mix(ALIGN_STAR.out.umi_dedup_tsv_per_umi)
        ch_umi_dedup_tsv_umi_per_position = ch_umi_dedup_tsv_umi_per_position.mix(ALIGN_STAR.out.umi_dedup_tsv_umi_per_position)
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.multiqc_files)
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    if (params.aligner == 'star_rsem') {

        QUANTIFY_RSEM (
            ch_transcriptome_bam,
            ch_rsem_index,
            params.use_sentieon_star
        )
        ch_rsem_stat              = QUANTIFY_RSEM.out.stat
        ch_rsem_logs              = QUANTIFY_RSEM.out.logs
        ch_rsem_counts_gene       = QUANTIFY_RSEM.out.counts_gene
        ch_rsem_counts_transcript = QUANTIFY_RSEM.out.counts_transcript
        ch_rsem_tpm_gene          = QUANTIFY_RSEM.out.merged_tpm_gene
        ch_rsem_tpm_transcript    = QUANTIFY_RSEM.out.merged_tpm_transcript
        ch_rsem_merged_counts_gene       = QUANTIFY_RSEM.out.merged_counts_gene
        ch_rsem_merged_counts_transcript = QUANTIFY_RSEM.out.merged_counts_transcript
        ch_rsem_merged_genes_long        = QUANTIFY_RSEM.out.merged_genes_long
        ch_rsem_merged_isoforms_long     = QUANTIFY_RSEM.out.merged_isoforms_long
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stat.collect{ tuple -> tuple[1] })
        ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_deseq2_pca   = ch_deseq2_pca.mix(DESEQ2_QC_RSEM.out.pca_multiqc)
            ch_deseq2_dists = ch_deseq2_dists.mix(DESEQ2_QC_RSEM.out.dists_multiqc)
            ch_deseq2_pdf          = ch_deseq2_pdf.mix(DESEQ2_QC_RSEM.out.pdf)
            ch_deseq2_rdata        = ch_deseq2_rdata.mix(DESEQ2_QC_RSEM.out.rdata)
            ch_deseq2_pca_txt      = ch_deseq2_pca_txt.mix(DESEQ2_QC_RSEM.out.pca_txt)
            ch_deseq2_dists_txt    = ch_deseq2_dists_txt.mix(DESEQ2_QC_RSEM.out.dists_txt)
            ch_deseq2_log          = ch_deseq2_log.mix(DESEQ2_QC_RSEM.out.log)
            ch_deseq2_size_factors = ch_deseq2_size_factors.mix(DESEQ2_QC_RSEM.out.size_factors)
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }

    } else if (params.aligner == 'star_salmon') {

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_STAR_SALMON (
            ch_samplesheet.map { item -> [ [:], item ] },
            ch_transcriptome_bam,
            ch_dummy_file,
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            true,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd
        )
        ch_star_salmon_quant         = QUANTIFY_STAR_SALMON.out.results
        ch_star_salmon_tx2gene       = QUANTIFY_STAR_SALMON.out.tx2gene
        ch_star_salmon_counts_gene   = QUANTIFY_STAR_SALMON.out.counts_gene
        ch_star_salmon_counts_gene_length_scaled = QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled
        ch_star_salmon_counts_gene_scaled = QUANTIFY_STAR_SALMON.out.counts_gene_scaled
        ch_star_salmon_counts_transcript = QUANTIFY_STAR_SALMON.out.counts_transcript
        ch_star_salmon_lengths_gene  = QUANTIFY_STAR_SALMON.out.lengths_gene
        ch_star_salmon_lengths_transcript = QUANTIFY_STAR_SALMON.out.lengths_transcript
        ch_star_salmon_tpm_gene      = QUANTIFY_STAR_SALMON.out.tpm_gene
        ch_star_salmon_tpm_transcript = QUANTIFY_STAR_SALMON.out.tpm_transcript
        ch_star_salmon_merged_gene_rds       = QUANTIFY_STAR_SALMON.out.merged_gene_rds_unified
        ch_star_salmon_merged_transcript_rds = QUANTIFY_STAR_SALMON.out.merged_transcript_rds_unified
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map { tuple -> tuple[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_deseq2_pca   = ch_deseq2_pca.mix(DESEQ2_QC_STAR_SALMON.out.pca_multiqc)
            ch_deseq2_dists = ch_deseq2_dists.mix(DESEQ2_QC_STAR_SALMON.out.dists_multiqc)
            ch_deseq2_pdf          = ch_deseq2_pdf.mix(DESEQ2_QC_STAR_SALMON.out.pdf)
            ch_deseq2_rdata        = ch_deseq2_rdata.mix(DESEQ2_QC_STAR_SALMON.out.rdata)
            ch_deseq2_pca_txt      = ch_deseq2_pca_txt.mix(DESEQ2_QC_STAR_SALMON.out.pca_txt)
            ch_deseq2_dists_txt    = ch_deseq2_dists_txt.mix(DESEQ2_QC_STAR_SALMON.out.dists_txt)
            ch_deseq2_log          = ch_deseq2_log.mix(DESEQ2_QC_STAR_SALMON.out.log)
            ch_deseq2_size_factors = ch_deseq2_size_factors.mix(DESEQ2_QC_STAR_SALMON.out.size_factors)
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { item -> [ [:], item ] },
            ch_splicesites.map { item -> [ [:], item ] },
            ch_fasta.map { item -> [ [:], item ] },
            params.with_umi,
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            params.bam_csi_index,
            params.skip_markduplicates,
            ch_transcriptome_bam,
            ch_transcript_fasta.map { item -> [ [:], item ] },
            ch_genome_bam,
            ch_genome_bam_index
        )
        ch_genome_bam          = ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = ALIGN_HISAT2.out.bai
        ch_unprocessed_bams    = ALIGN_HISAT2.out.orig_bam.map { meta, bam -> [ meta, bam, '' ] }
        ch_unaligned_sequences = ALIGN_HISAT2.out.unaligned
        ch_hisat2_summary      = ALIGN_HISAT2.out.summary
        ch_sorted_bam_stats    = ch_sorted_bam_stats.mix(ALIGN_HISAT2.out.stats)
        ch_sorted_bam_flagstat = ch_sorted_bam_flagstat.mix(ALIGN_HISAT2.out.flagstat)
        ch_sorted_bam_idxstats = ch_sorted_bam_idxstats.mix(ALIGN_HISAT2.out.idxstats)
        ch_umi_genomic_dedup_log        = ch_umi_genomic_dedup_log.mix(ALIGN_HISAT2.out.umi_genomic_dedup_log)
        ch_umi_transcriptomic_dedup_log = ch_umi_transcriptomic_dedup_log.mix(ALIGN_HISAT2.out.umi_transcriptomic_dedup_log)
        ch_umi_prepare_for_rsem_log     = ch_umi_prepare_for_rsem_log.mix(ALIGN_HISAT2.out.umi_prepare_for_rsem_log)
        ch_umi_transcriptome_dedup_bam      = ch_umi_transcriptome_dedup_bam.mix(ALIGN_HISAT2.out.umi_transcriptome_dedup_bam)
        ch_umi_transcriptome_sorted_bam     = ch_umi_transcriptome_sorted_bam.mix(ALIGN_HISAT2.out.umi_transcriptome_sorted_bam)
        ch_umi_transcriptome_sorted_bam_bai = ch_umi_transcriptome_sorted_bam_bai.mix(ALIGN_HISAT2.out.umi_transcriptome_sorted_bam_bai)
        ch_umi_transcriptome_filtered_bam   = ch_umi_transcriptome_filtered_bam.mix(ALIGN_HISAT2.out.umi_transcriptome_filtered_bam)
        ch_umi_dedup_stats    = ch_umi_dedup_stats.mix(ALIGN_HISAT2.out.umi_dedup_stats)
        ch_umi_dedup_bam      = ch_umi_dedup_bam.mix(ALIGN_HISAT2.out.umi_dedup_bam)
        ch_umi_dedup_bai      = ch_umi_dedup_bai.mix(ALIGN_HISAT2.out.umi_dedup_bai)
        ch_umi_dedup_flagstat = ch_umi_dedup_flagstat.mix(ALIGN_HISAT2.out.umi_dedup_flagstat)
        ch_umi_dedup_idxstats = ch_umi_dedup_idxstats.mix(ALIGN_HISAT2.out.umi_dedup_idxstats)
        ch_umi_dedup_tsv_edit_distance    = ch_umi_dedup_tsv_edit_distance.mix(ALIGN_HISAT2.out.umi_dedup_tsv_edit_distance)
        ch_umi_dedup_tsv_per_umi          = ch_umi_dedup_tsv_per_umi.mix(ALIGN_HISAT2.out.umi_dedup_tsv_per_umi)
        ch_umi_dedup_tsv_umi_per_position = ch_umi_dedup_tsv_umi_per_position.mix(ALIGN_HISAT2.out.umi_dedup_tsv_umi_per_position)
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_HISAT2.out.multiqc_files)
        ch_versions = ch_versions.mix(ALIGN_HISAT2.out.versions)
    }

    // Filter bam and index by percent mapped being present in the meta

    ch_genome_bam_bai_mapping = ch_genome_bam
        .join(ch_genome_bam_index)
        .join(ch_percent_mapped, remainder: true)
        .map{ row ->
            def (meta, bam, index) = row[0..2]
            def percent_mapped = row.size() == 4 ? row[3] : null
            def pass = percent_mapped != null ? percent_mapped >= params.min_mapped_reads.toFloat() : null
            return [ meta, bam, index, percent_mapped, pass ]
        }
        .multiMap { meta, bam, index, percent_mapped, pass ->
            bam: [ meta, bam, index, pass ]
            percent_mapped: [ meta.id, percent_mapped ]
            percent_mapped_pass: [ meta.id, percent_mapped, pass ]
            status: [ meta.id, pass ]
        }

    ch_percent_mapped = ch_genome_bam_bai_mapping.percent_mapped

    // Save mapping status for workflow summary where present

    ch_map_status = ch_genome_bam_bai_mapping.status
        .filter { _id, pass -> pass != null }

    // Save status for MultiQC report
    ch_fail_mapping_multiqc = ch_genome_bam_bai_mapping.percent_mapped_pass
        .filter { _id, _percent_mapped, pass -> pass != null && !pass }
        .map { id, percent_mapped, _pass -> [ "${id}\t${percent_mapped}" ] }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "STAR uniquely mapped reads (%)"]
                sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
        }
        .collectFile(name: 'fail_mapped_samples_mqc.tsv')

    ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc)

    // Where a percent mapping is present, use it to filter bam and index

    map_filtered_genome_bam_bai = ch_genome_bam_bai_mapping.bam
        .filter { _meta, _bam, _index, pass -> pass || pass == null }
        .multiMap { meta, bam, index, _pass ->
            bam: [ meta, bam ]
            index: [ meta, index ]
        }

    ch_genome_bam = map_filtered_genome_bam_bai.bam
    ch_genome_bam_index = map_filtered_genome_bam_bai.index

    //
    // MODULE: Run Preseq
    //
    if (!params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_txt = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_preseq_log = PRESEQ_LCEXTRAP.out.log
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{ tuple -> tuple[1] })
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    if (!params.skip_markduplicates && !params.with_umi) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            ch_fasta.map { item -> [ [:], item ] },
            ch_fai.map { item -> [ [:], item ] }
        )
        ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index = params.bam_csi_index ? BAM_MARKDUPLICATES_PICARD.out.csi : BAM_MARKDUPLICATES_PICARD.out.bai
        ch_markdup_bam      = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_markdup_bai      = params.bam_csi_index ? BAM_MARKDUPLICATES_PICARD.out.csi : BAM_MARKDUPLICATES_PICARD.out.bai
        ch_markdup_metrics  = BAM_MARKDUPLICATES_PICARD.out.metrics
        ch_markdup_stats    = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_markdup_flagstat = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_markdup_idxstats = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{ tuple -> tuple[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{ tuple -> tuple[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{ tuple -> tuple[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{ tuple -> tuple[1] })

        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam,
            ch_gtf
        )
        ch_stringtie_gtf        = STRINGTIE_STRINGTIE.out.transcript_gtf
        ch_stringtie_coverage   = STRINGTIE_STRINGTIE.out.coverage_gtf
        ch_stringtie_abundance  = STRINGTIE_STRINGTIE.out.abundance
        ch_stringtie_ballgown   = STRINGTIE_STRINGTIE.out.ballgown
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    if (!params.skip_qc && !params.skip_biotype_qc && biotype) {

        ch_gtf
            .map { gtf -> biotypeInGtf(gtf, biotype) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(ch_gtf)
            .combine(biotype_in_gtf)
            .filter { tuple -> tuple[-1] }
            .map { tuple -> tuple[0..<tuple.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_featurecounts_counts  = SUBREAD_FEATURECOUNTS.out.counts
        ch_featurecounts_summary = SUBREAD_FEATURECOUNTS.out.summary
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_biotype_counts = MULTIQC_CUSTOM_BIOTYPE.out.tsv
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{ tuple -> tuple[1] })
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    // Note: Strand parameters are conditional on library strandedness (see nextflow.config)
    //
    if (!params.skip_bigwig) {

        ch_genomecov_input = ch_genome_bam.map { meta, bam -> [ meta, bam, 1 ] }

        BEDTOOLS_GENOMECOV_FW (
            ch_genomecov_input,
            [],
            'bedGraph',
            true
        )
        BEDTOOLS_GENOMECOV_REV (
            ch_genomecov_input,
            [],
            'bedGraph',
            true
        )

        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_FW.out.versions.first())

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV_FW.out.genomecov,
            ch_chrom_sizes
        )
        ch_bigwig_forward = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.bigwig

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV_REV.out.genomecov,
            ch_chrom_sizes
        )
        ch_bigwig_reverse = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE.out.bigwig
    }

    //
    // MODULE: Downstream QC steps
    //
    if (!params.skip_qc) {
        if (!params.skip_qualimap) {
            // Sort BAM by name for qualimap (performance optimization)
            SAMTOOLS_SORT_QUALIMAP (
                ch_genome_bam,
                ch_fasta.map { item -> [ [:], item ] },
                ''
            )

            QUALIMAP_RNASEQ (
                SAMTOOLS_SORT_QUALIMAP.out.bam,
                ch_gtf.map { item -> [ [:], item ] }
            )
            ch_qualimap_results = QUALIMAP_RNASEQ.out.results
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.collect{ tuple -> tuple[1] })
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }

        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                ch_gtf.map { item -> [ [:], item ] }
            )
            ch_dupradar_scatter   = DUPRADAR.out.scatter2d
            ch_dupradar_boxplot   = DUPRADAR.out.boxplot
            ch_dupradar_histogram = DUPRADAR.out.hist
            ch_dupradar_gene_data = DUPRADAR.out.dupmatrix
            ch_dupradar_intercept = DUPRADAR.out.intercept_slope
            ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.collect{ tuple -> tuple[1] })
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        // Get RSeqC modules to run
        def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ module -> module.trim().toLowerCase() } : []
        if (params.bam_csi_index) {
            ['read_distribution', 'inner_distance', 'tin'].each { rseqc_module ->
                if (rseqc_modules.contains(rseqc_module)) {
                    rseqc_modules.remove(rseqc_module)
                }
            }
        }
        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]).map { meta, bam, bai -> [ meta, [ bam, bai ] ] },
                ch_gene_bed,
                rseqc_modules
            )
            ch_rseqc_bamstat              = BAM_RSEQC.out.bamstat_txt
            ch_rseqc_inferexperiment      = BAM_RSEQC.out.inferexperiment_txt
            ch_rseqc_junctionannotation_bed         = BAM_RSEQC.out.junctionannotation_bed
            ch_rseqc_junctionannotation_interact_bed = BAM_RSEQC.out.junctionannotation_interact_bed
            ch_rseqc_junctionannotation_xls          = BAM_RSEQC.out.junctionannotation_xls
            ch_rseqc_junctionannotation_log = BAM_RSEQC.out.junctionannotation_log
            ch_rseqc_junctionannotation_pdf        = BAM_RSEQC.out.junctionannotation_pdf
            ch_rseqc_junctionannotation_events_pdf = BAM_RSEQC.out.junctionannotation_events_pdf
            ch_rseqc_junctionannotation_r          = BAM_RSEQC.out.junctionannotation_rscript
            ch_rseqc_junctionsaturation_pdf = BAM_RSEQC.out.junctionsaturation_pdf
            ch_rseqc_junctionsaturation_r   = BAM_RSEQC.out.junctionsaturation_rscript
            ch_rseqc_readduplication_pos_xls = BAM_RSEQC.out.readduplication_pos_xls
            ch_rseqc_readduplication_seq_xls = BAM_RSEQC.out.readduplication_seq_xls
            ch_rseqc_readduplication_pdf  = BAM_RSEQC.out.readduplication_pdf
            ch_rseqc_readduplication_r    = BAM_RSEQC.out.readduplication_rscript
            ch_rseqc_readdistribution     = BAM_RSEQC.out.readdistribution_txt
            ch_rseqc_innerdistance_txt      = BAM_RSEQC.out.innerdistance_freq
            ch_rseqc_innerdistance_distance = BAM_RSEQC.out.innerdistance_distance
            ch_rseqc_innerdistance_mean     = BAM_RSEQC.out.innerdistance_mean
            ch_rseqc_innerdistance_pdf      = BAM_RSEQC.out.innerdistance_pdf
            ch_rseqc_innerdistance_r        = BAM_RSEQC.out.innerdistance_rscript
            ch_rseqc_tin                  = BAM_RSEQC.out.tin_txt
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.inferexperiment_txt.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionsaturation_rscript.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readduplication_pos_xls.collect{ tuple -> tuple[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.tin_txt.collect{ tuple -> tuple[1] })

            // Compare predicted supplied or Salmon-predicted strand with what we get from RSeQC
            ch_strand_comparison = BAM_RSEQC.out.inferexperiment_txt
                .map {
                    meta, strand_log ->
                        def rseqc_inferred_strand = getInferexperimentStrandedness(strand_log, params.stranded_threshold, params.unstranded_threshold)
                        def rseqc_strandedness = rseqc_inferred_strand.inferred_strandedness

                        def status = 'fail'
                        def multiqc_lines = []
                        if (meta.salmon_strand_analysis) {
                            def salmon_strandedness = meta.salmon_strand_analysis.inferred_strandedness

                            if (salmon_strandedness == rseqc_strandedness && rseqc_strandedness != 'undetermined') {
                                status = 'pass'
                            }
                            multiqc_lines = [
                                "$meta.id \tSalmon\t$status\tauto\t${meta.salmon_strand_analysis.values().join('\t')}",
                                "$meta.id\tRSeQC\t$status\tauto\t${rseqc_inferred_strand.values().join('\t')}"
                            ]
                        }
                        else {
                            if (meta.strandedness == rseqc_strandedness) {
                                status = 'pass'
                            }
                            multiqc_lines = [ "$meta.id\tRSeQC\t$status\t$meta.strandedness\t${rseqc_inferred_strand.values().join('\t')}" ]
                        }
                        return [ meta, status, multiqc_lines ]
                }
                .multiMap {
                    meta, status, multiqc_lines ->
                        status: [ meta.id, status == 'pass' ]
                        multiqc_lines: multiqc_lines
                }

            // Store the statuses for output
            ch_strand_status = ch_strand_comparison.status

            // Take the lines formatted for MultiQC and output
            ch_strand_comparison.multiqc_lines
                .flatten()
                .collect()
                .map {
                    tsv_data ->
                        def header = [
                            "Sample",
                            "Strand inference method",
                            "Status",
                            "Provided strandedness",
                            "Inferred strandedness",
                            "Sense (%)",
                            "Antisense (%)",
                            "Unstranded (%)"
                        ]
                        sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }

            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv'))
        }

        if (params.contaminant_screening in ['kraken2', 'kraken2_bracken'] ) {
            KRAKEN2 (
                ch_unaligned_sequences,
                params.kraken_db,
                params.save_kraken_assignments,
                params.save_kraken_unassigned
            )
            ch_kraken_reports = KRAKEN2.out.report
            ch_kraken_report  = KRAKEN2.out.report
            ch_versions = ch_versions.mix(KRAKEN2.out.versions)

            if (params.contaminant_screening == 'kraken2') {
                ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{ tuple -> tuple[1] })
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                BRACKEN (
                    ch_kraken_reports,
                    params.kraken_db
                )
                ch_bracken_txt = BRACKEN.out.txt
                ch_versions = ch_versions.mix(BRACKEN.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(BRACKEN.out.txt.collect{ tuple -> tuple[1] })
            }
        } else if (params.contaminant_screening == 'sylph') {
            def sylph_databases = params.sylph_db ? params.sylph_db.split(',').collect{ path -> file(path.trim()) } : []
            ch_sylph_databases = channel.value(sylph_databases)
            SYLPH_PROFILE (
                ch_unaligned_sequences,
                ch_sylph_databases
            )
            def sylph_profile_filtered = SYLPH_PROFILE.out.profile_out.filter{ tuple -> !tuple[1].isEmpty() }
            ch_sylph_profile = SYLPH_PROFILE.out.profile_out
            ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions)

            def sylph_taxonomies = params.sylph_taxonomy ? params.sylph_taxonomy.split(',').collect{ path -> file(path.trim()) } : []
            ch_sylph_taxonomies = channel.value(sylph_taxonomies)
            SYLPHTAX_TAXPROF (
                sylph_profile_filtered,
                ch_sylph_taxonomies
            )
            ch_sylphtax_output = SYLPHTAX_TAXPROF.out.taxprof_output
            ch_versions = ch_versions.mix(SYLPHTAX_TAXPROF.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SYLPHTAX_TAXPROF.out.taxprof_output.collect{ tuple -> tuple[1] })
        }
    }

    //
    // SUBWORKFLOW: Pseudoalignment and quantification with Salmon
    //
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {

        if (params.pseudo_aligner == 'salmon') {
            ch_pseudo_index = ch_salmon_index
        } else {
            ch_pseudo_index = ch_kallisto_index
        }

        QUANTIFY_PSEUDO_ALIGNMENT (
            ch_samplesheet.map { item -> [ [:], item ] },
            ch_strand_inferred_filtered_fastq,
            ch_pseudo_index,
            ch_dummy_file,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.pseudo_aligner,
            false,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd
        )
        ch_counts_gene_length_scaled = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled
        ch_pseudo_quant              = QUANTIFY_PSEUDO_ALIGNMENT.out.results
        ch_pseudo_tx2gene            = QUANTIFY_PSEUDO_ALIGNMENT.out.tx2gene
        ch_pseudo_counts_gene        = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene
        ch_pseudo_counts_gene_length_scaled = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled
        ch_pseudo_counts_gene_scaled = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_scaled
        ch_pseudo_counts_transcript  = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_transcript
        ch_pseudo_lengths_gene       = QUANTIFY_PSEUDO_ALIGNMENT.out.lengths_gene
        ch_pseudo_lengths_transcript = QUANTIFY_PSEUDO_ALIGNMENT.out.lengths_transcript
        ch_pseudo_tpm_gene           = QUANTIFY_PSEUDO_ALIGNMENT.out.tpm_gene
        ch_pseudo_tpm_transcript     = QUANTIFY_PSEUDO_ALIGNMENT.out.tpm_transcript
        ch_pseudo_merged_gene_rds       = QUANTIFY_PSEUDO_ALIGNMENT.out.merged_gene_rds_unified
        ch_pseudo_merged_transcript_rds = QUANTIFY_PSEUDO_ALIGNMENT.out.merged_transcript_rds_unified
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc.collect{ tuple -> tuple[1] })
        ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_PSEUDO (
                ch_counts_gene_length_scaled.map { tuple -> tuple[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            // Use separate channels for pseudo-aligner DESeq2 (published to different path)
            ch_pseudo_deseq2_pca   = ch_pseudo_deseq2_pca.mix(DESEQ2_QC_PSEUDO.out.pca_multiqc)
            ch_pseudo_deseq2_dists = ch_pseudo_deseq2_dists.mix(DESEQ2_QC_PSEUDO.out.dists_multiqc)
            ch_pseudo_deseq2_pdf          = ch_pseudo_deseq2_pdf.mix(DESEQ2_QC_PSEUDO.out.pdf)
            ch_pseudo_deseq2_rdata        = ch_pseudo_deseq2_rdata.mix(DESEQ2_QC_PSEUDO.out.rdata)
            ch_pseudo_deseq2_pca_txt      = ch_pseudo_deseq2_pca_txt.mix(DESEQ2_QC_PSEUDO.out.pca_txt)
            ch_pseudo_deseq2_dists_txt    = ch_pseudo_deseq2_dists_txt.mix(DESEQ2_QC_PSEUDO.out.dists_txt)
            ch_pseudo_deseq2_log          = ch_pseudo_deseq2_log.mix(DESEQ2_QC_PSEUDO.out.log)
            ch_pseudo_deseq2_size_factors = ch_pseudo_deseq2_size_factors.mix(DESEQ2_QC_PSEUDO.out.size_factors)
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_PSEUDO.out.versions)
        }
    }

    //
    // Collate and save software versions
    // Combines traditional versions.yml files with versions emitted via topic channels
    //
    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()
    ch_multiqc_data   = channel.empty()
    ch_multiqc_plots  = channel.empty()

    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config        = channel.fromPath("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config) : channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? channel.fromPath(params.multiqc_logo)   : channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = channel.value(
            paramsSummaryMultiqc(
                paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            )
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = channel.value(
            methodsDescriptionText(
                params.multiqc_methods_description
                    ? file(params.multiqc_methods_description)
                    : file("$projectDir/workflows/rnaseq/assets/multiqc/methods_description_template.yml", checkIfExists: true)
            )
        ).collectFile(name: 'methods_description_mqc.yaml')

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        // Provide MultiQC with rename patterns to ensure it uses sample names
        // for single-techrep samples not processed by CAT_FASTQ.
        //
        // We only add mappings when the FASTQ simpleName differs from the sample ID.
        // This prevents duplicate/conflicting mappings when multiple samples share
        // the same FASTQ filename in different directories (see #1657).
        //
        // Note: _raw/_trimmed suffixes are handled via extra_fn_clean_exts in multiqc_config.yml
        ch_name_replacements = ch_fastq
            .map{ meta, reads ->
                def paired = reads[0][1] as boolean
                def suffixes = paired ? ['_1', '_2'] : ['']
                def mappings = []

                def fastq1_simplename = file(reads[0][0]).simpleName
                if (fastq1_simplename != meta.id) {
                    mappings << [fastq1_simplename, "${meta.id}${suffixes[0]}"]
                    if (paired) {
                        mappings << [file(reads[0][1]).simpleName, "${meta.id}${suffixes[1]}"]
                    }
                }

                return mappings.collect { mapping -> mapping.join('\t') }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)
            .ifEmpty([])

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_multiqc_data   = MULTIQC.out.data
        ch_multiqc_plots  = MULTIQC.out.plots
    }

    //
    // Generate samplesheet with BAM paths for future runs
    //

    if (!params.skip_alignment && params.save_align_intermeds) {
        // Create channel with original input info and BAM paths
        ch_fastq.map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_unprocessed_bams.map { meta, genome_bam, transcriptome_bam -> [ meta.id, meta, genome_bam, transcriptome_bam ] })
            .join(ch_percent_mapped)
            .transpose()
            .map { _id, _fastq_meta, reads, meta, genome_bam, transcriptome_bam, percent_mapped ->

                // Handle BAM paths (same for all runs of this sample)
                def genome_bam_published = meta.has_genome_bam ?
                    (meta.original_genome_bam ?: '') :
                    mapBamToPublishedPath(genome_bam, meta.id, params.aligner)

                def transcriptome_bam_published = meta.has_transcriptome_bam ?
                    (meta.original_transcriptome_bam ?: '') :
                    mapBamToPublishedPath(transcriptome_bam, meta.id, params.aligner)

                def fastq_1 = reads[0].toUriString()
                def fastq_2 = reads.size() > 1 ? reads[1].toUriString() : ''
                def mapped = percent_mapped != null ? percent_mapped : ''

                return "${meta.id},${fastq_1},${fastq_2},${meta.strandedness},${genome_bam_published},${mapped},${transcriptome_bam_published}"
            }
            .collectFile(
                name: 'samplesheet_with_bams.csv',
                storeDir: "${workflow.outputDir}/samplesheets",
                newLine: true,
                seed: 'sample,fastq_1,fastq_2,strandedness,genome_bam,percent_mapped,transcriptome_bam'
            )
    }

    emit:
    trim_status       = ch_trim_status       // channel: [id, boolean]
    map_status        = ch_map_status        // channel: [id, boolean]
    strand_status     = ch_strand_status     // channel: [id, boolean]
    multiqc_report    = ch_multiqc_report    // channel: /path/to/multiqc_report.html
    multiqc_data      = ch_multiqc_data      // channel: /path/to/multiqc_data/
    multiqc_plots     = ch_multiqc_plots     // channel: /path/to/multiqc_plots/
    versions          = ch_versions          // channel: [ path(versions.yml) ]
    collated_versions = ch_collated_versions // channel: /path/to/collated_versions.yml

    // QC and trimming outputs
    fastqc_raw_html    = ch_fastqc_raw_html
    fastqc_raw_zip     = ch_fastqc_raw_zip
    fastqc_trim_html   = ch_fastqc_trim_html
    fastqc_trim_zip    = ch_fastqc_trim_zip
    trim_html          = ch_trim_html
    trim_zip           = ch_trim_zip
    trim_log           = ch_trim_log
    trim_json          = ch_trim_json
    trim_unpaired      = ch_trim_unpaired
    umi_log            = ch_umi_log
    umi_reads          = ch_umi_reads
    umi_genomic_dedup_log        = ch_umi_genomic_dedup_log
    umi_transcriptomic_dedup_log = ch_umi_transcriptomic_dedup_log
    umi_prepare_for_rsem_log     = ch_umi_prepare_for_rsem_log
    umi_transcriptome_dedup_bam      = ch_umi_transcriptome_dedup_bam
    umi_transcriptome_sorted_bam     = ch_umi_transcriptome_sorted_bam
    umi_transcriptome_sorted_bam_bai = ch_umi_transcriptome_sorted_bam_bai
    umi_transcriptome_filtered_bam   = ch_umi_transcriptome_filtered_bam
    umi_dedup_stats    = ch_umi_dedup_stats
    umi_dedup_bam      = ch_umi_dedup_bam
    umi_dedup_bai      = ch_umi_dedup_bai
    umi_dedup_flagstat = ch_umi_dedup_flagstat
    umi_dedup_idxstats = ch_umi_dedup_idxstats
    umi_dedup_tsv_edit_distance    = ch_umi_dedup_tsv_edit_distance
    umi_dedup_tsv_per_umi          = ch_umi_dedup_tsv_per_umi
    umi_dedup_tsv_umi_per_position = ch_umi_dedup_tsv_umi_per_position
    lint_log_raw       = ch_lint_log_raw
    lint_log_trimmed   = ch_lint_log_trimmed
    lint_log_bbsplit   = ch_lint_log_bbsplit
    lint_log_ribo      = ch_lint_log_ribo
    bbsplit_stats      = ch_bbsplit_stats
    sortmerna_log      = ch_sortmerna_log
    ribodetector_log   = ch_ribodetector_log
    seqkit_stats       = ch_seqkit_stats
    bowtie2_rrna_log   = ch_bowtie2_rrna_log
    bowtie2_rrna_index = ch_bowtie2_rrna_index
    seqkit_prefixed    = ch_seqkit_prefixed
    seqkit_converted   = ch_seqkit_converted

    // Alignment outputs
    star_log           = ch_star_log
    star_log_out       = ch_star_log_out
    star_log_progress  = ch_star_log_progress
    star_tab           = ch_star_tab
    star_bam           = ch_star_bam
    star_bai           = ch_star_bai
    sorted_bam_stats    = ch_sorted_bam_stats
    sorted_bam_flagstat = ch_sorted_bam_flagstat
    sorted_bam_idxstats = ch_sorted_bam_idxstats
    transcriptome_bam  = ch_transcriptome_bam_out
    unaligned_sequences = ch_unaligned_sequences
    hisat2_summary     = ch_hisat2_summary
    samtools_bai       = ch_samtools_bai  // Input BAM indices for BAM input mode

    // MarkDuplicates outputs
    markdup_bam        = ch_markdup_bam
    markdup_bai        = ch_markdup_bai
    markdup_metrics    = ch_markdup_metrics
    markdup_stats      = ch_markdup_stats
    markdup_flagstat   = ch_markdup_flagstat
    markdup_idxstats   = ch_markdup_idxstats

    // QC outputs
    preseq_txt         = ch_preseq_txt
    preseq_log         = ch_preseq_log
    qualimap_results   = ch_qualimap_results
    dupradar_scatter   = ch_dupradar_scatter
    dupradar_boxplot   = ch_dupradar_boxplot
    dupradar_histogram = ch_dupradar_histogram
    dupradar_gene_data = ch_dupradar_gene_data
    dupradar_intercept = ch_dupradar_intercept

    // RSeQC outputs
    rseqc_bamstat              = ch_rseqc_bamstat
    rseqc_inferexperiment      = ch_rseqc_inferexperiment
    rseqc_junctionannotation_bed         = ch_rseqc_junctionannotation_bed
    rseqc_junctionannotation_interact_bed = ch_rseqc_junctionannotation_interact_bed
    rseqc_junctionannotation_xls          = ch_rseqc_junctionannotation_xls
    rseqc_junctionannotation_log   = ch_rseqc_junctionannotation_log
    rseqc_junctionannotation_pdf        = ch_rseqc_junctionannotation_pdf
    rseqc_junctionannotation_events_pdf = ch_rseqc_junctionannotation_events_pdf
    rseqc_junctionannotation_r          = ch_rseqc_junctionannotation_r
    rseqc_junctionsaturation_pdf   = ch_rseqc_junctionsaturation_pdf
    rseqc_junctionsaturation_r     = ch_rseqc_junctionsaturation_r
    rseqc_readduplication_pos_xls  = ch_rseqc_readduplication_pos_xls
    rseqc_readduplication_seq_xls  = ch_rseqc_readduplication_seq_xls
    rseqc_readduplication_pdf      = ch_rseqc_readduplication_pdf
    rseqc_readduplication_r        = ch_rseqc_readduplication_r
    rseqc_readdistribution         = ch_rseqc_readdistribution
    rseqc_innerdistance_txt        = ch_rseqc_innerdistance_txt
    rseqc_innerdistance_distance   = ch_rseqc_innerdistance_distance
    rseqc_innerdistance_mean       = ch_rseqc_innerdistance_mean
    rseqc_innerdistance_pdf        = ch_rseqc_innerdistance_pdf
    rseqc_innerdistance_r          = ch_rseqc_innerdistance_r
    rseqc_tin                      = ch_rseqc_tin

    // Contaminant screening outputs
    kraken_report      = ch_kraken_report
    bracken_txt        = ch_bracken_txt
    sylph_profile      = ch_sylph_profile
    sylphtax_output    = ch_sylphtax_output

    // StringTie outputs → ${params.aligner}/stringtie
    stringtie_outputs = ch_stringtie_gtf
        .mix(ch_stringtie_coverage)
        .mix(ch_stringtie_abundance)
        .mix(ch_stringtie_ballgown)

    // FeatureCounts outputs → ${params.aligner}/featurecounts
    featurecounts_outputs = ch_featurecounts_counts
        .mix(ch_featurecounts_summary)
        .mix(ch_biotype_counts)

    // BigWig outputs → ${params.aligner}/bigwig
    bigwig_outputs = ch_bigwig_forward
        .mix(ch_bigwig_reverse)

    // Pseudo-alignment outputs → ${params.pseudo_aligner}
    pseudo_outputs = ch_pseudo_quant
        .mix(ch_pseudo_tx2gene)
        .mix(ch_pseudo_counts_gene)
        .mix(ch_pseudo_counts_gene_length_scaled)
        .mix(ch_pseudo_counts_gene_scaled)
        .mix(ch_pseudo_counts_transcript)
        .mix(ch_pseudo_lengths_gene)
        .mix(ch_pseudo_lengths_transcript)
        .mix(ch_pseudo_tpm_gene)
        .mix(ch_pseudo_tpm_transcript)
        .mix(ch_pseudo_merged_gene_rds)
        .mix(ch_pseudo_merged_transcript_rds)

    // RSEM outputs → star_rsem (logs to star_rsem/log)
    rsem_logs = ch_rsem_logs
    rsem_results = ch_rsem_stat
        .mix(ch_rsem_counts_gene)
        .mix(ch_rsem_counts_transcript)
        .mix(ch_rsem_tpm_gene)
        .mix(ch_rsem_tpm_transcript)
        .mix(ch_rsem_merged_counts_gene)
        .mix(ch_rsem_merged_counts_transcript)
        .mix(ch_rsem_merged_genes_long)
        .mix(ch_rsem_merged_isoforms_long)

    // STAR-Salmon outputs → star_salmon
    star_salmon_outputs = ch_star_salmon_quant
        .mix(ch_star_salmon_tx2gene)
        .mix(ch_star_salmon_counts_gene)
        .mix(ch_star_salmon_counts_gene_length_scaled)
        .mix(ch_star_salmon_counts_gene_scaled)
        .mix(ch_star_salmon_counts_transcript)
        .mix(ch_star_salmon_lengths_gene)
        .mix(ch_star_salmon_lengths_transcript)
        .mix(ch_star_salmon_tpm_gene)
        .mix(ch_star_salmon_tpm_transcript)
        .mix(ch_star_salmon_merged_gene_rds)
        .mix(ch_star_salmon_merged_transcript_rds)

    // DESeq2 outputs (aligner-based) → ${params.aligner}/deseq2_qc
    deseq2_outputs = ch_deseq2_pdf
        .mix(ch_deseq2_rdata)
        .mix(ch_deseq2_pca_txt)
        .mix(ch_deseq2_dists_txt)
        .mix(ch_deseq2_log)
        .mix(ch_deseq2_size_factors)

    // DESeq2 outputs (pseudo-aligner) → ${params.pseudo_aligner}/deseq2_qc
    pseudo_deseq2_outputs = ch_pseudo_deseq2_pdf
        .mix(ch_pseudo_deseq2_rdata)
        .mix(ch_pseudo_deseq2_pca_txt)
        .mix(ch_pseudo_deseq2_dists_txt)
        .mix(ch_pseudo_deseq2_log)
        .mix(ch_pseudo_deseq2_size_factors)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
