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
include { ALIGN_STAR                                        } from '../../subworkflows/local/align_star'
include { QUANTIFY_RSEM                                     } from '../../subworkflows/local/quantify_rsem'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_STAR_SALMON } from '../../subworkflows/local/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/local/quantify_pseudo_alignment'

include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList             } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getSalmonInferredStrandedness  } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getStarPercentMapped           } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { biotypeInGtf                   } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getInferexperimentStrandedness } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ               } from '../../modules/nf-core/cat/fastq'
include { BBMAP_BBSPLIT           } from '../../modules/nf-core/bbmap/bbsplit'
include { DUPRADAR                } from '../../modules/nf-core/dupradar'
include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort'
include { PRESEQ_LCEXTRAP         } from '../../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ         } from '../../modules/nf-core/qualimap/rnaseq'
include { STRINGTIE_STRINGTIE     } from '../../modules/nf-core/stringtie/stringtie'
include { SUBREAD_FEATURECOUNTS   } from '../../modules/nf-core/subread/featurecounts'
include { MULTIQC                 } from '../../modules/nf-core/multiqc'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../modules/nf-core/umitools/prepareforrsem'
include { SORTMERNA                                            } from '../../modules/nf-core/sortmerna'
include { SORTMERNA as SORTMERNA_INDEX                         } from '../../modules/nf-core/sortmerna'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-validation'
include { fromSamplesheet                  } from 'plugin/nf-validation'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../../subworkflows/nf-core/fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../../subworkflows/nf-core/fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_SORT_STATS_SAMTOOLS          } from '../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../../subworkflows/nf-core/bam_rseqc'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/workflows/rnaseq/assets/multiqc/biotypes_header.txt", checkIfExists: true)
ch_dummy_file                = ch_pca_header_multiqc

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
    ch_sortmerna_index   // channel: path(sortmerna/index/)
    ch_splicesites       // channel: path(genome.splicesites.txt)
    make_sortmerna_index // boolean: Whether to create an index before running sortmerna

    main:

    ch_multiqc_files = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            checkSamplesAfterGrouping(it)
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    ch_filtered_reads  = Channel.empty()
    ch_trim_read_count = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_filtered_reads  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_cat_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads  = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }
    ch_multiqc_files = ch_multiqc_files.mix(ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv'))

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            ch_bbsplit_index,
            [],
            [ [], [] ],
            false
        )
        .primary_fastq
        .set { ch_filtered_reads }
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    // Check rRNA databases for sortmerna
    if (params.remove_ribo_rna) {
        ch_ribo_db = file(params.ribo_database_manifest)
        if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}

        Channel.from(ch_ribo_db.readLines())
            .map { row -> file(row, checkIfExists: true) }
            .collect()
            .map { [ 'rrna_refs', it ] }
            .set { ch_sortmerna_fastas }

        if (make_sortmerna_index) {
            SORTMERNA_INDEX (
                [ [],[] ],
                ch_sortmerna_fastas,
                [ [],[] ]
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
        }

        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas,
            ch_sortmerna_index
        )
        .reads
        .set { ch_filtered_reads }

        ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log.collect{it[1]})
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_filtered_reads
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    ch_fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    def prepare_tool_indices = []
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {
        prepare_tool_indices << params.pseudo_aligner
    }
    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        !params.salmon_index && !('salmon' in prepare_tool_indices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_filtered_fastq }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam       = Channel.empty()
    ch_genome_bam_index = Channel.empty()
    ch_star_log         = Channel.empty()
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
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_star_log          = ALIGN_STAR.out.log_final
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
        ch_genome_bam       = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.bai
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
                    multiqcTsvFromList(tsv_data, header)
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

            BAM_RSEQC
                .out
                .inferexperiment_txt
                .map {
                    meta, strand_log ->
                        def inferred_strand = getInferexperimentStrandedness(strand_log, 30)
                        if (meta.strandedness != inferred_strand[0]) {
                            return [ "$meta.id\t$meta.strandedness\t${inferred_strand.join('\t')}" ]
                        }
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = [
                            "Sample",
                            "Provided strandedness",
                            "Inferred strandedness",
                            "Sense (%)",
                            "Antisense (%)",
                            "Undetermined (%)"
                        ]
                        multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }
            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv'))
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
        ch_multiqc_config        = Channel.fromPath("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
