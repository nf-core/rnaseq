/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowRnaseq.initialise(params, log)

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if (params.bam_csi_index) {
    for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
        if (rseqc_modules.contains(rseqc_module)) {
            rseqc_modules.remove(rseqc_module)
        }
    }
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../modules/local/deseq2_qc'
include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/local/umitools_prepareforrsem'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME                                    } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR                                        } from '../subworkflows/local/align_star'
include { QUANTIFY_RSEM                                     } from '../subworkflows/local/quantify_rsem'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_pseudo'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../subworkflows/local/quantify_pseudo'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit'
include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ             } from '../modules/nf-core/qualimap/rnaseq'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna'
include { STRINGTIE_STRINGTIE         } from '../modules/nf-core/stringtie/stringtie'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_HISAT2               } from '../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../subworkflows/nf-core/bam_rseqc'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow RNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
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
        is_aws_igenome,
        biotype,
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { WorkflowRnaseq.checkMaxContigSize(it, log) }
    }

    //
    // Create input channel from input file provided through params.input
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
            WorkflowRnaseq.validateInput(it)
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

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    PREPARE_GENOME.out.fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        !params.salmon_index && !('salmon' in prepareToolIndices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: WorkflowRnaseq.getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_strand_inferred_fastq,
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
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            PREPARE_GENOME.out.bbsplit_index,
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
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .set { ch_filtered_reads }

        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index.map { [ [:], it ] },
            PREPARE_GENOME.out.gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
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
            ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)

            // Co-ordinate sort, index and run stats on transcriptome BAM
            BAM_SORT_STATS_SAMTOOLS (
                ch_transcriptome_bam,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] }
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
                BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
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
                ch_umitools_dedup_bam.paired_end
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
            ch_transcriptome_bam,
            ch_dummy_file,
            PREPARE_GENOME.out.transcript_fasta,
            PREPARE_GENOME.out.gtf,
            'salmon',
            true,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_STAR_SALMON.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_STAR_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        QUANTIFY_RSEM (
            ch_filtered_reads,
            PREPARE_GENOME.out.rsem_index,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index  = QUANTIFY_RSEM.out.bai
        ch_samtools_stats    = QUANTIFY_RSEM.out.stats
        ch_samtools_flagstat = QUANTIFY_RSEM.out.flagstat
        ch_samtools_idxstats = QUANTIFY_RSEM.out.idxstats
        ch_star_multiqc      = QUANTIFY_RSEM.out.logs
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
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
            ch_aligner_pca_multiqc        = DESEQ2_QC_RSEM.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_RSEM.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_filtered_reads,
            PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] },
            PREPARE_GENOME.out.splicesites.map { [ [:], it ] },
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary
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
            ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
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
                    pass_mapped_reads[meta.id] = true
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    pass_mapped_reads[meta.id] = false
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
    }

    //
    // MODULE: Run Preseq
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_markduplicates && !params.with_umi) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fai.map { [ [:], it ] }
        )
        ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
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
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    ch_featurecounts_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

        PREPARE_GENOME
            .out
            .gtf
            .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(PREPARE_GENOME.out.gtf)
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
        ch_featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV.out.bedgraph_forward,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.versions)

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV.out.bedgraph_reverse,
            PREPARE_GENOME.out.chrom_sizes
        )
    }

    //
    // MODULE: Downstream QC steps
    //
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
    ch_tin_multiqc                = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
            QUALIMAP_RNASEQ (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf.map { [ [:], it ] }
            )
            ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }

        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_dupradar_multiqc = DUPRADAR.out.multiqc
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules
            )
            ch_bamstat_multiqc            = BAM_RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = BAM_RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = BAM_RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = BAM_RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = BAM_RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = BAM_RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = BAM_RSEQC.out.readduplication_pos_xls
            ch_tin_multiqc                = BAM_RSEQC.out.tin_txt
            ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

            ch_inferexperiment_multiqc
                .map {
                    meta, strand_log ->
                        def inferred_strand = WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30)
                        pass_strand_check[meta.id] = true
                        if (meta.strandedness != inferred_strand[0]) {
                            pass_strand_check[meta.id] = false
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
                        WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }
        }
    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_pseudo_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    
    if (!params.skip_pseudo_alignment) {

       if (params.pseudo_aligner == 'salmon'){
           ch_pseudo_index = PREPARE_GENOME.out.salmon_index
       } else {
           ch_pseudo_index = PREPARE_GENOME.out.kallisto_index
       }

        QUANTIFY_PSEUDO_ALIGNMENT (
            ch_filtered_reads,
            ch_pseudo_index,
            ch_dummy_file,
            PREPARE_GENOME.out.gtf,
            params.pseudo_aligner,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_pseudo_multiqc = QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc
        ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.versions)
        ch_counts_gene_length_scaled = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_PSEUDO (
                ch_counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_PSEUDO.out.pca_multiqc
            ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_PSEUDO.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_PSEUDO.out.versions)
        }
    }
    
    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
            ch_multiqc_logo.collect().ifEmpty([]),
            ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
            ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
            ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').ifEmpty([]),
            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_pseudo_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_aligner_pca_multiqc.collect().ifEmpty([]),
            ch_aligner_clustering_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_tin_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
    }
    
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)

    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
