/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
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
include { ALIGN_STAR    } from '../../subworkflows/local/align_star'
include { QUANTIFY_RSEM } from '../../subworkflows/local/quantify_rsem'
include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList             } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { getStarPercentMapped           } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { biotypeInGtf                   } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getInferexperimentStrandedness } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { methodsDescriptionText         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { DUPRADAR                   } from '../../modules/nf-core/dupradar'
include { SAMTOOLS_SORT              } from '../../modules/nf-core/samtools/sort'
include { PRESEQ_LCEXTRAP            } from '../../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ            } from '../../modules/nf-core/qualimap/rnaseq'
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie'
include { SUBREAD_FEATURECOUNTS      } from '../../modules/nf-core/subread/featurecounts'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { MULTIQC                    } from '../../modules/nf-core/multiqc'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../modules/nf-core/umitools/prepareforrsem'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_SORT_STATS_SAMTOOLS          } from '../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../../subworkflows/nf-core/bam_rseqc'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
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

// Header files for MultiQC
ch_pca_header_multiqc           = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
sample_status_header_multiqc    = file("$projectDir/workflows/rnaseq/assets/multiqc/sample_status_header.txt", checkIfExists: true)
ch_clustering_header_multiqc    = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc      = file("$projectDir/workflows/rnaseq/assets/multiqc/biotypes_header.txt", checkIfExists: true)
ch_dummy_file                   = ch_pca_header_multiqc

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
    ch_splicesites       // channel: path(genome.splicesites.txt)
    make_sortmerna_index // boolean: Whether to create an index before running sortmerna

    main:

    ch_multiqc_files = Channel.empty()
    ch_trim_status = Channel.empty()
    ch_map_status = Channel.empty()
    ch_strand_status = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            checkSamplesAfterGrouping(samplesheet)
        }
        .set { ch_fastq }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !salmon_index_available,
        !params.sortmerna_index && params.remove_ribo_rna,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        ALIGN_STAR (
            ch_strand_inferred_filtered_fastq,
            ch_star_index.map { [ [:], it ] },
            ch_gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = ALIGN_STAR.out.bam
        ch_genome_bam_index    = ALIGN_STAR.out.bai
        ch_transcriptome_bam   = ALIGN_STAR.out.bam_transcript
        ch_star_log            = ALIGN_STAR.out.log_final
        ch_unaligned_sequences = ALIGN_STAR.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            // Deduplicate genome BAM file before downstream analysis
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_genome_bam       = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.deduplog.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats.collect{it[1]})

            if (params.bam_csi_index) {
                ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)

            // Co-ordinate sort, index and run stats on transcriptome BAM
            BAM_SORT_STATS_SAMTOOLS (
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )
            ch_transcriptome_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_transcriptome_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

            // Deduplicate transcriptome BAM file before read counting with Salmon
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
                params.umitools_dedup_stats
            )

            // Name sort BAM before passing to Salmon
            SAMTOOLS_SORT (
                BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam,
                ch_fasta.map { [ [:], it ] }
            )

            // Only run prepare_for_rsem.py on paired-end BAM files
            SAMTOOLS_SORT
                .out
                .bam
                .branch {
                    meta, bam ->
                        single_end: meta.single_end
                            return [ meta, bam ]
                        paired_end: !meta.single_end
                            return [ meta, bam ]
                }
                .set { ch_umitools_dedup_bam }

            // Fix paired-end reads in name sorted BAM file
            // See: https://github.com/nf-core/rnaseq/issues/828
            UMITOOLS_PREPAREFORSALMON (
                ch_umitools_dedup_bam.paired_end.map { meta, bam -> [ meta, bam, [] ] }
            )
            ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

            ch_umitools_dedup_bam
                .single_end
                .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
                .set { ch_transcriptome_bam }
        }

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_STAR_SALMON (
            ch_samplesheet.map { [ [:], it ] },
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
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map { it[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        QUANTIFY_RSEM (
            ch_strand_inferred_filtered_fastq,
            ch_rsem_index,
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam       = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index = QUANTIFY_RSEM.out.bai
        ch_star_log         = QUANTIFY_RSEM.out.logs
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stat.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = QUANTIFY_RSEM.out.csi
        }
        ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { [ [:], it ] },
            ch_splicesites.map { [ [:], it ] },
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = FASTQ_ALIGN_HISAT2.out.bai
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_genome_bam       = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.deduplog.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats.collect{it[1]})
            if (params.bam_csi_index) {
                ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_log
            .map { meta, align_log -> [ meta ] + getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        // Save status for workflow summary
        ch_map_status = ch_percent_mapped
            .map {
                meta, mapped, pass ->
                    return [ meta.id, pass ]
            }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv'))
    }

    //
    // MODULE: Run Preseq
    //
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{it[1]})
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    if (!params.skip_alignment && !params.skip_markduplicates && !params.with_umi) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            ch_fasta.map { [ [:], it ] },
            ch_fai.map { [ [:], it ] }
        )
        ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam,
            ch_gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

        ch_gtf
            .map { biotypeInGtf(it, biotype) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(ch_gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{it[1]})
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

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
        ch_versions = ch_versions.mix(BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.versions)

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV_REV.out.genomecov,
            ch_chrom_sizes
        )
    }

    //
    // MODULE: Downstream QC steps
    //
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
            QUALIMAP_RNASEQ (
                ch_genome_bam,
                ch_gtf.map { [ [:], it ] }
            )
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.collect{it[1]})
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }

        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                ch_gtf.map { [ [:], it ] }
            )
            ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.collect{it[1]})
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        // Get RSeqC modules to run
        def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        if (params.bam_csi_index) {
            for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
                if (rseqc_modules.contains(rseqc_module)) {
                    rseqc_modules.remove(rseqc_module)
                }
            }
        }
        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_gene_bed,
                rseqc_modules
            )
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.inferexperiment_txt.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionsaturation_rscript.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readduplication_pos_xls.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.tin_txt.collect{it[1]})
            ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

            // Compare predicted supplied or Salmon-predicted strand with what we get from RSeQC
            ch_strand_comparison = BAM_RSEQC.out.inferexperiment_txt
                .map {
                    meta, strand_log ->
                        def rseqc_inferred_strand = getInferexperimentStrandedness(strand_log, params.stranded_threshold, params.unstranded_threshold)
                        rseqc_strandedness = rseqc_inferred_strand.inferred_strandedness

                        def status = 'fail'
                        def multiqc_lines = []
                        if (meta.salmon_strand_analysis) {
                            salmon_strandedness = meta.salmon_strand_analysis.inferred_strandedness

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
            ch_versions = ch_versions.mix(KRAKEN2.out.versions)

            if (params.contaminant_screening == 'kraken2') {
                ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{it[1]})
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                BRACKEN (
                    ch_kraken_reports,
                    params.kraken_db
                )
                ch_versions = ch_versions.mix(BRACKEN.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(BRACKEN.out.txt.collect{it[1]})
            }
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
            ch_samplesheet.map { [ [:], it ] },
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
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc.collect{it[1]})
        ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_PSEUDO (
                ch_counts_gene_length_scaled.map { it[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_PSEUDO.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()

    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config        = Channel.fromPath("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = Channel.value(
            paramsSummaryMultiqc(
                paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            )
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = Channel.value(
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
        // for single-techrep samples not processed by CAT_FASTQ, and trims out
        // _raw or _trimmed

        ch_name_replacements = ch_fastq
            .map{ meta, reads ->
                def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
                def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
                if (reads[0][1] ){
                    def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
                    def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
                    def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
                    return [ name1, name2, fastqcnames1, fastqcnames2 ]
                } else{
                    return [ name1, fastqcnames ]
                }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    trim_status    = ch_trim_status    // channel: [id, boolean]
    map_status     = ch_map_status     // channel: [id, boolean]
    strand_status  = ch_strand_status  // channel: [id, boolean]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
