/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC                          } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../../modules/local/deseq2_qc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../../modules/local/multiqc_custom_biotype'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { STAR_SALMON                    } from './star_salmon.nf'
include { STAR_RSEM                      } from './star_rsem.nf'

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
include { BEDTOOLS_GENOMECOV                                   } from '../../modules/nf-core/bedtools/genomecov'
include { UCSC_BEDCLIP                                         } from '../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG                                } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-validation'
include { fromSamplesheet                  } from 'plugin/nf-validation'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../../subworkflows/nf-core/fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIM } from './trim.nf'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../../subworkflows/nf-core/fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_SORT_STATS_SAMTOOLS          } from '../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../../subworkflows/nf-core/bam_rseqc'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/nf-core/quantify_pseudo_alignment'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    input                // : List<FastqRun>
    ch_fasta             // file: genome.fasta
    ch_gtf               // file: genome.gtf
    ch_fai               // file: genome.fai
    ch_chrom_sizes       // file: genome.sizes
    ch_gene_bed          // file: gene.bed
    ch_transcript_fasta  // file: transcript.fasta
    ch_star_index        // directory: star/index/
    ch_rsem_index        // directory: rsem/index/
    ch_hisat2_index      // directory: hisat2/index/
    ch_salmon_index      // directory: salmon/index/
    ch_kallisto_index    // directory: kallisto/index/
    ch_bbsplit_index     // directory: bbsplit/index/
    ch_sortmerna_index   // directory: sortmerna/index/
    ch_splicesites       // file: genome.splicesites.txt
    make_sortmerna_index // boolean: Whether to create an index before running sortmerna

    main:

    // Header files for MultiQC
    ch_pca_header_multiqc        = file("$moduleDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
    ch_clustering_header_multiqc = file("$moduleDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
    ch_biotypes_header_multiqc   = file("$moduleDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)
    ch_dummy_file                = ch_pca_header_multiqc

    //
    // Create channel from input file provided through params.input
    //
    ch_fastq = input
        .collect { run ->
            !run.fastq_2
                ? [ [ id:run.id, single_end:true ], [ run.fastq_1 ] ]
                : [ [ id:run.id, single_end:false ], [ run.fastq_1, run.fastq_2 ] ]
        }
        .groupBy { meta, fastqs -> meta.id }
        .values()
        .collect { runs ->
            checkSamplesAfterGrouping(runs)
            return [ runs.first().v1, runs.collectMany { _meta, fqs -> fqs } ]
        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    ch_cat_fastq = ch_fastq.map { meta, fastqs ->
        fastqs.size() > 1
            ? [ meta, CAT_FASTQ( meta, fastqs ) ]
            : [ meta, fastqs ]
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    ch_trim_reads = FASTQ_FASTQC_UMITOOLS_TRIM (
        ch_cat_fastq,
        params.trimmer,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_umi_extract,
        params.umi_discard_read,
        params.skip_trimming,
        null,
        params.save_trimmed_fail,
        params.save_merged
    )

    ch_filtered_reads = ch_trim_reads
        .filter { meta, fastqc, umi, trim ->
            trim.num_reads == null || trim.num_reads >= params.min_trimmed_reads.toLong()
        }
        .map { meta, fastqc, umi, trim -> [ meta, trim.reads ] }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_fail_trimming_records = ch_trim_reads
        .map { meta, fastqc, umi, trim -> [ meta, trim.num_reads ] }
        .filter { meta, num_reads -> num_reads < params.min_trimmed_reads.toLong() }
        .map { meta, num_reads -> [ meta.id, num_reads ] }
        .collect()
    ch_fail_trimming_multiqc = file('fail_trimmed_samples_mqc.tsv')
    mergeCsv(ch_fail_trimming_multiqc, ch_fail_trimming_records, header: ['Sample', 'Reads after trimming'], sep: '\t')

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        ch_filtered_reads = ch_filtered_reads.map { meta, fastq ->
            def bbsplit = BBMAP_BBSPLIT ( meta, fastq, ch_bbsplit_index, null, [], [], false )
            [ meta, bbsplit.primary_fastq ]
        }
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    // Check rRNA databases for sortmerna
    if (params.remove_ribo_rna) {
        ch_ribo_db = file(params.ribo_database_manifest)
        if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}

        ch_sortmerna_fastas = ch_ribo_db.readLines()
            .collect { row -> file(row, checkIfExists: true) }

        if (make_sortmerna_index) {
            ch_sortmerna_index = SORTMERNA_INDEX ( [:], [], ch_sortmerna_fastas, null ).index
        }

        ch_filtered_reads = ch_filtered_reads.map { meta, fastq ->
            def sortmerna = SORTMERNA ( meta, fastq, ch_sortmerna_fastas, ch_sortmerna_index )
            [ meta, sortmerna.reads ]
        }
    }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //

    def prepare_tool_indices = []
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {
        prepare_tool_indices << params.pseudo_aligner
    }

    ch_strand_inferred_filtered_fastq = FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_filtered_reads,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        !params.salmon_index && !('salmon' in prepare_tool_indices)
    )

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_bam = Channel.empty()
    ch_counts = null
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        (ch_bam, ch_counts) = STAR_SALMON (
            ch_strand_inferred_filtered_fastq,
            ch_star_index,
            ch_gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_fasta,
            ch_transcript_fasta,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.salmon_quant_libtype ?: '',
        )
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        (ch_bam, ch_counts) = STAR_RSEM (
            ch_strand_inferred_filtered_fastq,
            ch_rsem_index,
            ch_fasta
        )
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        (ch_bam, ch_counts) = FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta
        )
    }

    //
    // MODULE: Differential expression analysis with DeSeq2
    //
    if (!params.skip_qc & !params.skip_deseq2_qc) {
        // TODO: choice type for counts
        // star_salmon -> ch_counts.counts_gene_length_scaled
        // star_rsem -> ch_counts.merged_counts_gene
        deseq2 = DESEQ2_QC (
            ch_counts.gene,
            ch_pca_header_multiqc,
            ch_clustering_header_multiqc
        )
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_percent_mapped = ch_bam
            .map { sample -> [ sample, getStarPercentMapped(params, sample.align_log) ] }

        ch_bam = ch_percent_mapped
            .filter { sample, percent, pass -> pass }
            .map { sample, percent, pass -> sample }

        def records = ch_percent_mapped
            .filterMap { sample, percent, pass ->
                !pass ? [ sample.id, percent ] : null
            }
            .collect()

        ch_fail_mapping_multiqc = file('fail_mapped_samples_mqc.tsv')
        mergeCsv(ch_fail_mapping_multiqc, records, header: ["Sample", "STAR uniquely mapped reads (%)"], sep: '\t')
    }

    //
    // MODULE: Run Preseq
    //
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        preseq = ch_bam.map { sample ->
            PRESEQ_LCEXTRAP ( sample )
        }
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    if (!params.skip_alignment && !params.skip_markduplicates && !params.with_umi) {
        picard = BAM_MARKDUPLICATES_PICARD (
            ch_bam,
            ch_fasta,
            ch_fai
        )
        ch_genome_bam       = picard.markdup.filterMap { meta, out -> out.bam ? [ meta, out.bam ] : null }
        ch_genome_bam_index = picard.indexed.filterMap { meta, out ->
            def idx = params.bam_csi_index ? out.csi : out.bai
            idx ? [ meta, idx ] : null
        }
    }

    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
        stringtie = ch_genome_bam.map { meta, bam ->
            def args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
            def out = STRINGTIE_STRINGTIE ( meta, bam, ch_gtf, args )
            [ meta, out ]
        }
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        featurecounts = ch_genome_bam
            .filter { meta, bam -> biotypeInGtf(ch_gtf, biotype) }
            .map { meta, bam ->
                def args = [
                    '-B -C',
                    params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type",
                    "-t $params.featurecounts_feature_type"
                ].join(' ').trim()

                def (counts, summary) = SUBREAD_FEATURECOUNTS ( meta, [ bam ], ch_gtf, args )
                def biotype = MULTIQC_CUSTOM_BIOTYPE ( counts, ch_biotypes_header_multiqc )
                [ meta, counts, summary, biotype ]
            }
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

        bigwig_forward = ch_genome_bam.map { meta, bam ->
            def prefix = meta.strandedness == 'reverse' ? "${meta.id}.reverse" : "${meta.id}.forward"
            def args = '-split -du -strand + -bg'
            def genomecov = BEDTOOLS_GENOMECOV ( meta, bam, 1, null, 'bedGraph', true, args, prefix )
            def clipped = UCSC_BEDCLIP ( meta, clipped, ch_chrom_sizes, '', "${meta.id}.clip.forward" )
            def bigwig = UCSC_BEDGRAPHTOBIGWIG ( meta, clipped, ch_chrom_sizes, '', "${meta.id}.forward" )
            [ meta, bigwig ]
        }

        bigwig_reverse = ch_genome_bam.map { meta, bam ->
            def prefix = meta.strandedness == 'reverse' ? "${meta.id}.forward" : "${meta.id}.reverse"
            def args = '-split -du -strand - -bg'
            def genomecov = BEDTOOLS_GENOMECOV ( meta, bam, 1, null, 'bedGraph', true, args, prefix )
            def clipped = UCSC_BEDCLIP ( meta, clipped, ch_chrom_sizes, '', "${meta.id}.clip.reverse" )
            def bigwig = UCSC_BEDGRAPHTOBIGWIG ( meta, clipped, ch_chrom_sizes, '', "${meta.id}.reverse" )
            [ meta, bigwig ]
        }
    }

    //
    // MODULE: Downstream QC steps
    //
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
            qualimap = ch_genome_bam.map { meta, bam ->
                [ meta, QUALIMAP_RNASEQ ( meta, bam, ch_gtf ) ]
            }
        }

        if (!params.skip_dupradar) {
            dupradar = ch_genome_bam.map { meta, bam ->
                [ meta, DUPRADAR ( meta, bam, ch_gtf ) ]
            }
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
            rseqc = BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_gene_bed,
                rseqc_modules
            )

            def records = rseqc.inferexperiment
                .filterMap { meta, strand_log ->
                    def inferred_strand = getInferexperimentStrandedness(strand_log, 30)
                    if (meta.strandedness != inferred_strand[0]) {
                        return [ meta.id, meta.strandedness ] + inferred_strand
                    }
                    else {
                        return null
                    }
                }
                .collect()

            def header = [
                "Sample",
                "Provided strandedness",
                "Inferred strandedness",
                "Sense (%)",
                "Antisense (%)",
                "Undetermined (%)"
            ]
            ch_fail_strand_multiqc = file('fail_strand_check_mqc.tsv')
            mergeCsv(ch_fail_strand_multiqc, records, header: header, sep: '\t')
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

        quant = QUANTIFY_PSEUDO_ALIGNMENT (
            samplesheet,
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

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_PSEUDO (
                quant.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
        }
    }

    //
    // Collate and save software versions
    //
    collated_versions = file("${params.outdir}/pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    mergeText(collated_versions, softwareVersionsToYAML(), newLine: true)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config        = file("$moduleDir/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : null
        ch_multiqc_logo          = params.multiqc_logo   ? file(params.multiqc_logo)   : null
        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        workflow_summary         = paramsSummaryMultiqc(summary_params)
        mergeText('workflow_summary_mqc.yaml', workflow_summary)

        multiqc = MULTIQC (
            Channel.topic('logs').collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config,
            ch_multiqc_logo,
            params.multiqc_title ? "--title \"$params.multiqc_title\"" : ''
        )
        ch_multiqc_report = multiqc.report
    }

    emit:
    samples = ch_bam
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html

    // publish:
    // ch_bam          >> 'aligned'
    // ch_counts       >> 'counts'
    // deseq2          >> 'deseq2_qc'
    // preseq          >> 'preseq'
    // stringtie       >> 'stringtie'
    // featurecounts   >> 'featurecounts'
    // bigwig_forward  >> 'bigwig'
    // bigwig_reverse  >> 'bigwig'
    // qualimap        >> 'qualimap'
    // dupradar        >> 'dupradar'
    // picard          >> 'aligner'
    // rseqc           >> 'rseqc'
}

record FastqRun {
    id: String
    fastq_1: Path
    fastq_2: Path?
    strandedness: Strandedness
}

enum Strandedness {
    FORWARD,
    REVERSE,
    UNSTRANDED,
    AUTO
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
