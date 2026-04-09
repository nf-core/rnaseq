/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC as DESEQ2_QC_BAM_SALMON } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../../modules/local/deseq2_qc'
include { RUSTQC                              } from '../../modules/nf-core/rustqc/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ALIGN_STAR                            } from '../../subworkflows/local/align_star'
include { ALIGN_BOWTIE2                         } from '../../subworkflows/local/align_bowtie2'
include { BAM_POST_ALIGNMENT_QC                 } from '../../subworkflows/local/bam_post_alignment_qc'
include { QUANTIFY_RSEM                         } from '../../subworkflows/nf-core/quantify_rsem'
include { BAM_DEDUP_UMI                         } from '../../subworkflows/nf-core/bam_dedup_umi'

include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList             } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
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
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { SYLPH_PROFILE              } from '../../modules/nf-core/sylph/profile/main'
include { SYLPHTAX_TAXPROF           } from '../../modules/nf-core/sylphtax/taxprof/main'
include { MULTIQC                    } from '../../modules/nf-core/multiqc'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_BAM_SALMON } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet          // channel: path(sample_sheet.csv)
    ch_versions             // channel: [ path(versions.yml) ]
    ch_fasta                // channel: path(genome.fasta)
    ch_gtf                  // channel: path(genome.gtf)
    ch_fai                  // channel: path(genome.fai)
    ch_chrom_sizes          // channel: path(genome.sizes)
    ch_gene_bed             // channel: path(gene.bed)
    ch_transcript_fasta     // channel: path(transcript.fasta)
    ch_star_index           // channel: path(star/index/)
    ch_rsem_index           // channel: path(rsem/index/)
    ch_hisat2_index         // channel: path(hisat2/index/)
    ch_bowtie2_index        // channel: path(bowtie2/index/) for alignment
    ch_salmon_index         // channel: path(salmon/index/)
    ch_kallisto_index       // channel: [ meta, path(kallisto/index/) ]
    ch_bbsplit_index        // channel: path(bbsplit/index/)
    ch_ribo_db              // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index      // channel: path(sortmerna/index/)
    ch_bowtie2_rrna_index   // channel: path(bowtie2/index/) for rRNA removal
    ch_splicesites          // channel: path(genome.splicesites.txt)
    ch_kraken_db            // channel: path(kraken2/db/)

    main:

    // Header files for MultiQC
    def ch_pca_header_multiqc        = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
    def sample_status_header_multiqc = file("$projectDir/workflows/rnaseq/assets/multiqc/sample_status_header.txt", checkIfExists: true)
    def ch_clustering_header_multiqc = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
    def ch_biotypes_header_multiqc   = file("$projectDir/workflows/rnaseq/assets/multiqc/biotypes_header.txt", checkIfExists: true)
    def ch_transcript_fasta_placeholder = ch_pca_header_multiqc

    // Pre-build fasta_fai value channels for subworkflows that need [meta, fasta, fai]
    // .first() converts the queue channel to a value channel so it can be consumed multiple times
    ch_fasta_fai            = ch_fasta.combine(ch_fai).map { fasta, fai -> [ [:], fasta, fai ] }.first()
    ch_transcript_fasta_fai = ch_transcript_fasta.map { fasta -> [[:], fasta, []] }

    ch_multiqc_files = channel.empty()
    ch_trim_status = channel.empty()
    ch_map_status = channel.empty()
    ch_strand_status = channel.empty()
    ch_percent_mapped = channel.empty()
    ch_unaligned_sequences = channel.empty()

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
    ch_genome_bam_index = SAMTOOLS_INDEX.out.index

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
        ch_bowtie2_rrna_index,                      // ch_bowtie2_index (for rRNA removal)
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
        params.save_merged_fastq,                   // save_merged_fastq
        params.stranded_threshold,                  // stranded_threshold
        params.unstranded_threshold                 // unstranded_threshold
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
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
    ch_star_log            = channel.empty()

    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem')) {
        // Use iGenomes-pinned STAR (2.6.1d) only when aligning with a pre-built iGenomes index
        def use_igenomes_star = params.genome && params.star_index ? true : false

        ALIGN_STAR (
            ch_strand_inferred_filtered_fastq,
            ch_star_index.map { item -> [ [:], item ] },
            ch_gtf.map { item -> [ [:], item ] },
            params.star_ignore_sjdbgtf,
            use_igenomes_star,
            ch_fasta_fai,
            params.use_sentieon_star,
            params.use_parabricks_star,
            params.skip_markduplicates
        )

        ch_genome_bam                    = ch_genome_bam.mix(ALIGN_STAR.out.bam)
        ch_genome_bam_index              = ch_genome_bam_index.mix(ALIGN_STAR.out.index)
        ch_transcriptome_bam             = ch_transcriptome_bam.mix(ALIGN_STAR.out.bam_transcript)
        ch_percent_mapped                = ch_percent_mapped.mix(ALIGN_STAR.out.percent_mapped)
        ch_unprocessed_bams              = ch_genome_bam.join(ch_transcriptome_bam)
        ch_star_log                      = ALIGN_STAR.out.log_final
        ch_unaligned_sequences           = ALIGN_STAR.out.fastq
        ch_multiqc_files                 = ch_multiqc_files.mix(ch_star_log)

        if (!params.with_umi && (params.skip_markduplicates || params.use_parabricks_star)) {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise. If mark duplicates
            // will run, those stats will be added later instead to avoid
            // duplicate flagstat files in MultiQC.
            // When Parabricks handles markduplicates internally, Picard is
            // skipped, so we also need to add alignment stats here.

            ch_multiqc_files = ch_multiqc_files
                .mix(ALIGN_STAR.out.stats)
                .mix(ALIGN_STAR.out.flagstat)
                .mix(ALIGN_STAR.out.idxstats)
        }
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2
    //
    ch_bowtie2_log = channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2_salmon') {

        ALIGN_BOWTIE2 (
            ch_strand_inferred_filtered_fastq,
            ch_bowtie2_index,
            ch_fasta_fai
        )

        // For Bowtie2+Salmon, the BAM is aligned to transcriptome so it's the "transcriptome_bam"
        // Use orig_bam (query-grouped) for Salmon - coordinate-sorted BAM breaks paired-end quantification
        ch_genome_bam                    = ch_genome_bam.mix(ALIGN_BOWTIE2.out.bam)
        ch_genome_bam_index              = ch_genome_bam_index.mix(ALIGN_BOWTIE2.out.index)
        ch_transcriptome_bam             = ch_transcriptome_bam.mix(ALIGN_BOWTIE2.out.orig_bam)
        ch_percent_mapped                = ch_percent_mapped.mix(ALIGN_BOWTIE2.out.percent_mapped)
        ch_unprocessed_bams              = ch_genome_bam.map { meta, bam -> [ meta, bam, '' ] }
        ch_bowtie2_log                   = ALIGN_BOWTIE2.out.log_final
        ch_multiqc_files                 = ch_multiqc_files.mix(ch_bowtie2_log)

        if (!params.with_umi && params.skip_markduplicates) {
            ch_multiqc_files = ch_multiqc_files
                .mix(ALIGN_BOWTIE2.out.stats)
                .mix(ALIGN_BOWTIE2.out.flagstat)
                .mix(ALIGN_BOWTIE2.out.idxstats)
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { item -> [ [:], item ] },
            ch_splicesites.map { item -> [ [:], item ] },
            ch_fasta_fai,
            params.save_unaligned || (params.contaminant_screening && params.contaminant_screening_input == 'unmapped')
        )
        ch_genome_bam          = ch_genome_bam.mix(FASTQ_ALIGN_HISAT2.out.bam)
        ch_genome_bam_index    = ch_genome_bam_index.mix(FASTQ_ALIGN_HISAT2.out.index)
        ch_unprocessed_bams    = ch_genome_bam.map { meta, bam -> [ meta, bam, '' ] }
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary)

        if (!params.with_umi && params.skip_markduplicates) {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise. If mark duplicates
            // will run, those stats will be added later instead to avoid
            // duplicate flagstat files in MultiQC.
            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats)
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat)
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats)
        }
    }

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //
    if (!params.skip_alignment && params.with_umi) {

        BAM_DEDUP_UMI(
            ch_genome_bam.join(ch_genome_bam_index, by: [0]),
            ch_fasta_fai,
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            ch_transcriptome_bam,
            ch_transcript_fasta_fai,
            params.umitools_dedup_primary_only
        )

        ch_genome_bam        = BAM_DEDUP_UMI.out.bam
        ch_transcriptome_bam = BAM_DEDUP_UMI.out.transcriptome_bam
        ch_genome_bam_index  = BAM_DEDUP_UMI.out.index

        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_DEDUP_UMI.out.multiqc_files)
    }

    //
    // Quantification
    //
    if (params.aligner == 'star_rsem') {

        QUANTIFY_RSEM (
            ch_samplesheet.map { item -> [ [:], item ] },
            ch_transcriptome_bam,
            ch_rsem_index,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.use_sentieon_star,
            params.skip_quantification_merge
        )
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stat)

        if (!params.skip_qc & !params.skip_deseq2_qc & !params.skip_quantification_merge) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.counts_gene_length_scaled.map { _meta, counts -> counts },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect().map { file -> [[:], file] })
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect().map { file -> [[:], file] })
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }

    } else if (params.aligner in ['star_salmon', 'bowtie2_salmon']) {

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_BAM_SALMON (
            ch_samplesheet.map { item -> [ [:], item ] },
            ch_transcriptome_bam,
            ch_transcript_fasta_placeholder,
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            true,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge
        )
        if (!params.skip_qc & !params.skip_deseq2_qc & !params.skip_quantification_merge) {
            DESEQ2_QC_BAM_SALMON (
                QUANTIFY_BAM_SALMON.out.counts_gene_length_scaled.map { _meta, counts -> counts },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_BAM_SALMON.out.pca_multiqc.collect().map { file -> [[:], file] })
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_BAM_SALMON.out.dists_multiqc.collect().map { file -> [[:], file] })
            ch_versions = ch_versions.mix(DESEQ2_QC_BAM_SALMON.out.versions)
        }
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

    ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.map { file -> [[:], file] })

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
    // SUBWORKFLOW: Mark duplicate reads
    //

    // Some tools (Ex. Parabricks) may have already run marked duplicates during alignment
    def markdups_done = !params.skip_markduplicates && params.use_parabricks_star
    if (!params.skip_markduplicates && !params.with_umi && !markdups_done) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            ch_fasta_fai
        )
        ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.index
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics)
    }

    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam,
            ch_gtf
        )
    }

    //
    // Pre-compute param-derived values for QC subworkflow
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ module -> module.trim().toLowerCase() } : []
    if (params.bam_csi_index) {
        rseqc_modules.removeAll(['read_distribution', 'inner_distance', 'tin'])
    }

    ch_inferexperiment_txt = channel.empty()

    if (!params.skip_qc) {
        if (params.use_rustqc) {
            //
            // MODULE: RustQC - single-pass replacement for multiple QC tools
            //
            RUSTQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_gtf.map { gtf -> [ [:], gtf ] },
            )

            // Collect all per-tool outputs, flatten to individual files, then filter for MultiQC
            ch_multiqc_files = ch_multiqc_files.mix(
                RUSTQC.out.dupradar
                    .mix(RUSTQC.out.featurecounts)
                    .mix(RUSTQC.out.preseq)
                    .mix(RUSTQC.out.samtools)
                    .mix(RUSTQC.out.rseqc)
                    .mix(RUSTQC.out.qualimap)
                    .flatMap { meta, files -> (files instanceof List ? files : [files]).collect { f -> [meta, f] } }
                    .filter { _meta, f ->
                        // Exclude gene-level featureCounts summary so MultiQC only sees the
                        // biotype-level summary (*.biotype.tsv.summary), matching the default
                        // pipeline's featureCounts -g gene_biotype output.
                        if (f.name.endsWith('.featureCounts.tsv.summary')) return false
                        f.name =~ /(?i)\.(txt|tsv|xls|log|stats|flagstat|idxstats|html)$/ || f.name.contains('_mqc.')
                    }
            )

            // Extract infer_experiment from rseqc channel
            ch_inferexperiment_txt = RUSTQC.out.rseqc
                .map { meta, files ->
                    def ie = (files instanceof List ? files : [files]).find { it.name.endsWith('.infer_experiment.txt') }
                    ie ? [meta, ie] : null
                }
                .filter { it != null }
        } else {
            //
            // SUBWORKFLOW: Post-alignment QC (upstream tools)
            //
            BAM_POST_ALIGNMENT_QC (
                ch_genome_bam,
                ch_genome_bam_index,
                ch_gtf,
                ch_gene_bed,
                ch_fasta_fai,
                ch_biotypes_header_multiqc,
                params.skip_preseq,
                params.skip_biotype_qc,
                params.skip_qualimap,
                params.skip_dupradar,
                params.skip_rseqc,
                biotype,
                rseqc_modules
            )
            ch_multiqc_files = ch_multiqc_files.mix(BAM_POST_ALIGNMENT_QC.out.multiqc_files)
            ch_inferexperiment_txt = BAM_POST_ALIGNMENT_QC.out.inferexperiment_txt
            ch_versions = ch_versions.mix(BAM_POST_ALIGNMENT_QC.out.versions)
        }

        //
        // Strandedness comparison using infer_experiment output
        //
        // RustQC always produces infer_experiment output regardless of rseqc_modules
        def run_infer_experiment = params.use_rustqc || rseqc_modules.contains('infer_experiment')
        if (run_infer_experiment) {

            // Compare predicted supplied or Salmon-predicted strand with what we get from RSeQC
            ch_strand_comparison = ch_inferexperiment_txt
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

            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').map { file -> [[:], file] })
        }

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

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV_FW.out.genomecov,
            ch_chrom_sizes
        )

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV_REV.out.genomecov,
            ch_chrom_sizes
        )
    }

    if (!params.skip_qc) {
        //
        // Contaminant screening (Kraken2/Bracken/Sylph)
        //
        def ch_contaminant_sequences = params.contaminant_screening_input == 'trimmed'
            ? ch_strand_inferred_filtered_fastq
            : ch_unaligned_sequences

        if (params.contaminant_screening in ['kraken2', 'kraken2_bracken'] ) {
            KRAKEN2 (
                ch_contaminant_sequences,
                ch_kraken_db,
                params.save_kraken_assignments,
                params.save_kraken_unassigned
            )
            ch_kraken_reports = KRAKEN2.out.report

            if (params.contaminant_screening == 'kraken2') {
                ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report)
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                BRACKEN (
                    ch_kraken_reports,
                    ch_kraken_db
                )
                ch_multiqc_files = ch_multiqc_files.mix(BRACKEN.out.txt)
            }
        } else if (params.contaminant_screening == 'sylph') {
            def sylph_databases = params.sylph_db ? params.sylph_db.split(',').collect{ path -> file(path.trim()) } : []
            ch_sylph_databases = channel.value(sylph_databases)
            SYLPH_PROFILE (
                ch_contaminant_sequences,
                ch_sylph_databases
            )
            ch_sylph_profile = SYLPH_PROFILE.out.profile_out.filter{ tuple -> !tuple[1].isEmpty() }

            def sylph_taxonomies = params.sylph_taxonomy ? params.sylph_taxonomy.split(',').collect{ path -> file(path.trim()) } : []
            ch_sylph_taxonomies = channel.value(sylph_taxonomies)
            SYLPHTAX_TAXPROF (
                ch_sylph_profile,
                ch_sylph_taxonomies
            )
            ch_multiqc_files = ch_multiqc_files.mix(SYLPHTAX_TAXPROF.out.taxprof_output)
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
            ch_transcript_fasta_placeholder,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.pseudo_aligner,
            false,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge
        )
        ch_counts_gene_length_scaled = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc)

        if (!params.skip_qc & !params.skip_deseq2_qc & !params.skip_quantification_merge) {
            DESEQ2_QC_PSEUDO (
                ch_counts_gene_length_scaled.map { _meta, counts -> counts },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.pca_multiqc.collect().map { file -> [[:], file] })
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.dists_multiqc.collect().map { file -> [[:], file] })
            ch_versions = ch_versions.mix(DESEQ2_QC_PSEUDO.out.versions)
        }
    }

    //
    // Collate and save software versions
    // Combines traditional versions.yml files with versions emitted via topic channels
    //
    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()

    if (!params.skip_multiqc) {

        // Prepare MultiQC config paths
        def mqc_default_config = file("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        def mqc_custom_config  = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []
        def mqc_logo           = params.multiqc_logo   ? file(params.multiqc_logo, checkIfExists: true)   : []

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

        // Add global context files (summary, versions, methods) with empty meta
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary.map   { file -> [[:], file] })
            .mix(ch_collated_versions.map  { file -> [[:], file] })
            .mix(ch_methods_description.map{ file -> [[:], file] })

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

        // Build the MultiQC config list (default + optional custom)
        def mqc_config_files = mqc_custom_config ? [mqc_default_config, mqc_custom_config] : [mqc_default_config]

        if (params.skip_quantification_merge) {
            //
            // Per-sample MultiQC: group files by sample, attach global files to each
            //
            ch_multiqc_branched = ch_multiqc_files
                .branch { meta, _file ->
                    per_sample: meta.id != null
                    global: true
                }

            ch_global_files = ch_multiqc_branched.global
                .map { _meta, file -> file }
                .collect()

            ch_multiqc_input = ch_multiqc_branched.per_sample
                .map { meta, file -> [meta.id, file] }
                .groupTuple()
                .combine(ch_global_files.toList())
                .map { id, sample_files, global_files ->
                    [
                        [id: id],
                        sample_files + (global_files ?: []),
                        mqc_config_files,
                        mqc_logo,
                        [],  // replace_names not needed: each per-sample report contains only one sample's files, so there is no filename ambiguity to resolve
                        [],  // no sample_names
                    ]
                }
        } else {
            //
            // Normal merged MultiQC: collect all files into one report
            // Note: meta.id = 'multiqc_report' is a sentinel used by
            // modules/nf-core/multiqc/nextflow.config to select the merged
            // output path and prefix
            //
            ch_all_mqc_files = ch_multiqc_files
                .map { _meta, file -> file }
                .collect()

            ch_multiqc_input = ch_all_mqc_files
                .combine(ch_name_replacements.ifEmpty([]).toList())
                .map { items ->
                    def files = items[0..-2]
                    def replace_names = items[-1] instanceof List ? [] : items[-1]
                    [
                        [id: 'multiqc_report'],
                        files,
                        mqc_config_files,
                        mqc_logo,
                        replace_names,
                        [],
                    ]
                }
        }

        MULTIQC(ch_multiqc_input)
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> report }
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
                    mapBamToPublishedPath(genome_bam, meta.id, params.aligner, params.outdir)

                def transcriptome_bam_published = meta.has_transcriptome_bam ?
                    (meta.original_transcriptome_bam ?: '') :
                    mapBamToPublishedPath(transcriptome_bam, meta.id, params.aligner, params.outdir)

                def fastq_1 = reads[0].toUriString()
                def fastq_2 = reads.size() > 1 ? reads[1].toUriString() : ''
                def mapped = percent_mapped != null ? percent_mapped : ''

                def seq_platform = meta.seq_platform ?: params.seq_platform ?: ''
                def seq_center = meta.seq_center ?: params.seq_center ?: ''

                return "${meta.id},${fastq_1},${fastq_2},${meta.strandedness},${seq_platform},${seq_center},${genome_bam_published},${mapped},${transcriptome_bam_published}"
            }
            .collectFile(
                name: 'samplesheet_with_bams.csv',
                storeDir: "${params.outdir}/samplesheets",
                newLine: true,
                seed: 'sample,fastq_1,fastq_2,strandedness,seq_platform,seq_center,genome_bam,percent_mapped,transcriptome_bam'
            )
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
