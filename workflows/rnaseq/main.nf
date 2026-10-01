nextflow.enable.types = true

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
include { MULTIQC_RNASEQ                        } from '../../subworkflows/local/multiqc_rnaseq'
include { BAM_QC_RNASEQ                         } from '../../subworkflows/nf-core/bam_qc_rnaseq'
include { QUANTIFY_RSEM                         } from '../../subworkflows/nf-core/quantify_rsem'
include { BAM_DEDUP_UMI                         } from '../../subworkflows/nf-core/bam_dedup_umi'

include { ReadsInput; SampleRow; StringtieInput } from '../../modules/nf-core/types'
include { Bowtie2Aligned; StarAligned; MultiqcFiles; AlignedSample; Bam; Contaminants; StringtieSample; BigwigSample; PipelineInfo; UmiDedupBam; MarkdupBam; BamQcRnaseq; StringtieMerged; Hisat2Aligned; RrnaReferences; FastqQcTrimFilterSetstrandedness; QuantMerged; SampleRuns; TrimReadCount; TrimStatus; PercentMapped; MapStatus; PercentMappedPass; InferExperimentLog; StrandData; StrandStatus } from '../../modules/nf-core/types'
include { RsemMergeSample } from '../../modules/nf-core/custom/rsemmergecounts/main'
include { KallistoQuantSample } from '../../modules/nf-core/kallisto/quant/main'
include { MultiqcReport } from '../../modules/nf-core/multiqc/main'
include { RsemQuantSample } from '../../modules/nf-core/rsem/calculateexpression/main'
include { RustqcResult } from '../../modules/nf-core/rustqc/main'
include { SalmonQuantSample } from '../../modules/nf-core/salmon/quant/main'
include { SamtoolsIndexResult } from '../../modules/nf-core/samtools/index/main'
include { StringtieResult } from '../../modules/nf-core/stringtie/stringtie/main'
include { Deseq2Qc } from '../../modules/local/deseq2_qc/main'

include { readSamplesheet; samplesheetRowsToCsv } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline/samplesheet'
include { classifyStrand                 } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getHisat2PercentMapped         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

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
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_COMBINED    } from '../../modules/nf-core/bedtools/genomecov'
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_STRINGTIE_MERGE              } from '../../subworkflows/nf-core/bam_stringtie_merge/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_COMBINED } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
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
    ch_sample_rows: Channel<SampleRow>                 // one row per sequencing run of each sample
    ch_fasta: Value<Path?>                             // genome.fasta
    ch_fai: Value<Path?>                               // genome.fai
    ch_gtf: Value<Path?>                               // genome.gtf
    ch_chrom_sizes: Value<Path?>                       // genome.sizes
    ch_gene_bed: Value<Path?>                          // gene.bed
    ch_transcript_fasta: Value<Path?>                  // transcript.fasta
    ch_star_index: Value<Path?>                        // star/index/
    ch_rsem_index: Value<Path?>                        // rsem/index/
    ch_hisat2_index: Value<Path?>                      // hisat2/index/
    ch_bowtie2_index: Value<Path?>                     // bowtie2/index/ for alignment
    ch_salmon_index: Value<Path?>                      // salmon/index/
    ch_kallisto_index: Value<Path?>                    // kallisto/index/
    ch_bbsplit_index: Value<Path?>                     // bbsplit/index/
    ch_ribo_db: Channel<Path>                          // sortmerna_fasta_list
    ch_sortmerna_index: Value<Path?>                   // sortmerna/index/
    ch_bowtie2_rrna_index: Value<Path?>                // bowtie2/index/ for rRNA removal
    ch_splicesites: Value<Path?>                       // genome.splicesites.txt
    ch_kraken_db: Value<Path?>                         // kraken2/db/
    qc_tools: List<String>                             // QC tools to run, e.g. ['preseq', 'qualimap', 'rseqc_bam_stat', ...]

    main:

    def ch_pca_header_multiqc        = file("$projectDir/assets/deseq2_pca_header.txt", checkIfExists: true)
    def sample_status_header_multiqc = file("$projectDir/assets/sample_status_header.txt", checkIfExists: true)
    def ch_clustering_header_multiqc = file("$projectDir/assets/deseq2_clustering_header.txt", checkIfExists: true)
    def ch_biotypes_header_multiqc   = file("$projectDir/assets/biotypes_header.txt", checkIfExists: true)

    // Match the General Statistics column the active aligner emits so the
    // MultiQC fail_mapped row reads consistently with the rest of the report.
    def aligner_display_name = [
        'star_salmon'    : 'STAR uniquely mapped reads',
        'star_rsem'      : 'STAR uniquely mapped reads',
        'hisat2'         : 'HISAT2 overall alignment rate',
        'bowtie2_salmon' : 'Bowtie2 overall alignment rate',
    ][params.aligner as String] ?: 'Aligned reads'

    // Files each stage contributes to MultiQC, per sample. The report-only and sample-only
    // channels hold the files that go to just the merged report or just the per-sample reports
    // (--skip_quantification_merge). fail_* rows are appended inside MULTIQC_RNASEQ.
    def ch_mqc_files: Channel<MultiqcFiles>       = channel.empty()
    def ch_mqc_sample_only: Channel<MultiqcFiles> = channel.empty()
    def ch_mqc_report_only: Channel<Path>         = channel.empty()

    // The rows are validated and merged per sample as one batch, since a sample's runs can be
    // spread over any rows.
    def ch_input = ch_sample_rows
        .collect()
        .flatMap { rows -> readSamplesheet(rows, "${projectDir}/assets/schema_input.json", params.skip_alignment) as List<Map> }
        .map { s ->
            record(
                id:                s.id,
                meta:              s.meta,
                reads:             s.reads,
                runs:              s.runs,
                bam:               s.bam,
                transcriptome_bam: s.transcriptome_bam,
                percent_mapped:    s.percent_mapped,
                prealigned:        s.prealigned
            )
        }

    // Samplesheet re-assembled from the validated rows, for the SummarizedExperiment sample metadata
    def ch_samplesheet: Value<Path> = ch_sample_rows
        .collect()
        .flatMap { rows -> [ samplesheetRowsToCsv(rows) ] }
        .collectFile(name: 'samplesheet.csv')
        .collect()
        .map { files -> files.toList().first() as Path }

    // Samples that go through FASTQ preprocessing and samples supplied as pre-aligned BAM files
    ch_fastq_samples = ch_input.filter { s -> !s.prealigned }
    ch_bam_samples   = ch_input.filter { s -> s.prealigned }

    // One entry per sequencing run of each sample
    def ch_fastq: Channel<SampleRuns> = ch_fastq_samples.map { s -> record(id: s.id, meta: s.meta, runs: s.runs) }

    // Index pre-aligned genome BAM files; a sample may supply only a transcriptome BAM
    ch_prealigned_genome = ch_bam_samples.filter { s -> s.bam != null }
    def ch_bam_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_prealigned_genome)
    ch_prealigned_indexed = ch_prealigned_genome.join(ch_bam_index, by: 'id', remainder: true)
    ch_prealigned_indexed.subscribe { r ->
        if( r.bam == null || r.bai == null ) {
            error "Sample '${r.id}' is missing its pre-aligned BAM index result"
        }
    }
    ch_prealigned = ch_prealigned_indexed.filter { r -> r.bam != null && r.bai != null }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The bowtie2 rRNA index is built inside the FASTQ subworkflow, not in PREPARE_GENOME_INDICES.
    def make_bowtie2_index = !params.bowtie2_rrna_index && params.remove_ribo_rna && params.ribo_removal_tool == 'bowtie2'

    fastq_preprocessed = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq_samples,                           // ch_reads
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
        false,                                      // make_salmon_index (PREPARE_GENOME_INDICES already builds/loads this)
        false,                                      // make_sortmerna_index (PREPARE_GENOME_INDICES already builds/loads this)
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

    def ch_preprocessed: Channel<FastqQcTrimFilterSetstrandedness> = fastq_preprocessed.samples

    // Run-level rRNA references, built by FASTQ_REMOVE_RRNA from the rRNA
    // FASTAs only when no bowtie2 rRNA index was supplied
    def ch_rrna_references: Value<RrnaReferences> = fastq_preprocessed.rrna_references

    // Samples that fail min_trimmed_reads have no filtered reads and go no further
    ch_reads_ok = ch_preprocessed.filter { r -> r.reads != null }

    // MultiQC files from FASTQ preprocessing. The trimmer log only reaches the merged
    // report when the trimmer is not fastp.
    ch_fastq_qc = ch_preprocessed.map { r ->
        def fastqc_head = (r.fastqc?.raw_zip ?: []) + (r.fastqc?.trim_zip ?: [])
        def trim_json   = r.trim?.json ?: []
        def trim_log    = r.trim?.log ?: []
        def other       = [r.umi?.log, r.bbsplit?.stats, r.rrna?.sortmerna_log, r.rrna?.ribodetector_log, r.rrna?.seqkit_stats, r.rrna?.bowtie2_log] +
                          (r.fastqc?.filtered_zip ?: [])
        record(
            id:    r.id,
            files: (fastqc_head + (params.trimmer == 'fastp' ? [] : trim_log) + trim_json + other).findAll { f -> f != null }
        )
    }
    ch_mqc_files = ch_mqc_files.mix(ch_fastq_qc)
    if (params.trimmer == 'fastp') {
        ch_mqc_sample_only = ch_mqc_sample_only.mix(
            ch_preprocessed.map { r -> record(id: r.id, files: r.trim?.log ?: []) }
        )
    }

    def ch_trim_read_count: Channel<TrimReadCount> = ch_preprocessed
        .filter { r -> r.num_trimmed_reads != null }
        .map { r -> record(id: r.id, meta: r.meta, num_reads: r.num_trimmed_reads) }

    def ch_trim_status: Channel<TrimStatus> = ch_preprocessed
        .filter { r -> r.num_trimmed_reads != null }
        .map { r -> record(id: r.id, pass: r.num_trimmed_reads > params.min_trimmed_reads.toFloat()) }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    def ch_star: Channel<StarAligned> = channel.empty()
    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem')) {
        ch_star = ALIGN_STAR (
            ch_reads_ok,
            ch_star_index,
            ch_gtf,
            params.star_ignore_sjdbgtf,
            ch_fasta,
            ch_fai,
            params.use_sentieon_star,
            params.use_parabricks_star,
            params.skip_markduplicates
        )

        ch_mqc_files = ch_mqc_files.mix(ch_star.map { r -> record(id: r.id, files: [r.star.log_final]) })

        if (!params.with_umi && (params.skip_markduplicates || params.use_parabricks_star)) {
            // When Picard markduplicates runs, its stats are added later; adding these too would
            // duplicate flagstat files in MultiQC. Parabricks skips Picard, so add them here.
            ch_mqc_files = ch_mqc_files.mix(
                ch_star.map { r -> record(id: r.id, files: [r.samtools.stats, r.samtools.flagstat, r.samtools.idxstats]) }
            )
        }
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2
    //
    def ch_bowtie2: Channel<Bowtie2Aligned> = channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2_salmon') {
        ch_bowtie2 = ALIGN_BOWTIE2 (
            ch_reads_ok,
            ch_bowtie2_index,
            ch_fasta,
            ch_fai
        )

        ch_mqc_files = ch_mqc_files.mix(ch_bowtie2.map { r -> record(id: r.id, files: [r.bowtie2.log]) })

        if (!params.with_umi && params.skip_markduplicates) {
            ch_mqc_files = ch_mqc_files.mix(
                ch_bowtie2.map { r -> record(id: r.id, files: [r.samtools.stats, r.samtools.flagstat, r.samtools.idxstats]) }
            )
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    def ch_hisat2: Channel<Hisat2Aligned> = channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        ch_hisat2 = FASTQ_ALIGN_HISAT2 (
            ch_reads_ok,
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta,
            ch_fai,
            params.save_unaligned || (params.contaminant_screening && params.contaminant_screening_input == 'unmapped')
        )

        ch_mqc_files = ch_mqc_files.mix(ch_hisat2.map { r -> record(id: r.id, files: [r.hisat2.summary]) })

        if (!params.with_umi && params.skip_markduplicates) {
            ch_mqc_files = ch_mqc_files.mix(
                ch_hisat2.map { r -> record(id: r.id, files: [r.samtools.stats, r.samtools.flagstat, r.samtools.idxstats]) }
            )
        }
    }

    // The aligner outputs have different shapes, so the channel is left untyped
    def ch_aligned = channel.empty().mix(ch_star).mix(ch_bowtie2).mix(ch_hisat2)

    //
    // Genome-aligned BAMs with their index and mapping percentage, from the aligner or the samplesheet
    //
    def ch_genome_bam: Channel<Bam> = ch_prealigned
        .mix(ch_star)
        .mix(ch_bowtie2)
        .mix(ch_hisat2.map { r -> r + record(percent_mapped: getHisat2PercentMapped(r.hisat2.summary)) })

    def ch_transcriptome_bam: Channel<Bam> = ch_star
        .mix(ch_bowtie2)
        .mix(ch_bam_samples)
        .filter { r -> r.transcriptome_bam != null }

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //
    def ch_umi_dedup: Channel<UmiDedupBam> = channel.empty()
    if (!params.skip_alignment && params.with_umi) {
        ch_umi_dedup = BAM_DEDUP_UMI(
            ch_genome_bam,
            ch_fasta,
            ch_fai,
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            ch_transcriptome_bam,
            ch_transcript_fasta,
            params.umitools_dedup_primary_only
        )

        // The right-hand record wins on every shared field, so the deduplicated bam, bai and samtools replace the aligner's
        ch_genome_deduped = ch_genome_bam.join(ch_umi_dedup, by: 'id', remainder: true)
        ch_genome_deduped.subscribe { r ->
            if( r.genomic_dedup_log == null ) {
                error "Sample '${r.id}' is missing its UMI deduplication result"
            }
        }
        ch_genome_bam = ch_genome_deduped.filter { r -> r.genomic_dedup_log != null }
        ch_transcriptome_bam = ch_umi_dedup.filter { r -> r.transcriptome_bam != null }

        // Genome-side files only; MultiQC cannot tell transcriptome stats apart from genome stats
        ch_mqc_files = ch_mqc_files.mix(
            ch_umi_dedup.map { r -> record(id: r.id, files: [r.genomic_dedup_log, r.samtools.stats, r.samtools.flagstat, r.samtools.idxstats]) }
        )
    }

    //
    // Quantification
    //
    def ch_transcriptome_reads: Channel<ReadsInput> = ch_transcriptome_bam.map { r -> record(id: r.id, meta: r.meta, reads: [r.transcriptome_bam]) }

    def run_deseq2_qc = !params.skip_qc && !params.skip_deseq2_qc && !params.skip_quantification_merge

    def ch_quant: Channel<RsemQuantSample>             = channel.empty()
    def ch_quant_salmon: Channel<SalmonQuantSample>    = channel.empty()
    def ch_quant_merged: Channel<QuantMerged>          = channel.empty()
    def ch_quant_rsem_merge: Channel<RsemMergeSample>  = channel.empty()
    def ch_deseq2: Channel<Deseq2Qc>                   = channel.empty()
    if (params.aligner == 'star_rsem') {
        rsem = QUANTIFY_RSEM (
            ch_samplesheet,
            ch_transcriptome_reads,
            ch_rsem_index,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.use_sentieon_star,
            params.skip_quantification_merge
        )

        ch_quant            = rsem.samples
        ch_quant_merged     = rsem.merged
        ch_quant_rsem_merge = rsem.rsem_merge

        ch_mqc_files = ch_mqc_files.mix(rsem.samples.map { r -> record(id: r.id, files: [r.stat]) })

        if (run_deseq2_qc) {
            ch_deseq2 = DESEQ2_QC_RSEM (
                ch_quant_merged,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
        }
    } else if (params.aligner in ['star_salmon', 'bowtie2_salmon']) {

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        bam_salmon = QUANTIFY_BAM_SALMON (
            ch_samplesheet,
            ch_transcriptome_reads,
            null,
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge
        )

        ch_quant_salmon = bam_salmon.salmon
        ch_quant_merged = bam_salmon.merged

        if (run_deseq2_qc) {
            ch_deseq2 = DESEQ2_QC_BAM_SALMON (
                ch_quant_merged,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
        }
    }

    ch_mapped = ch_genome_bam.map { r ->
        r + record(pass: r.percent_mapped != null ? r.percent_mapped >= params.min_mapped_reads.toFloat() : null)
    }

    def ch_percent_mapped: Channel<PercentMapped> = ch_mapped.map { r -> record(id: r.id, percent_mapped: r.percent_mapped) }

    def ch_map_status: Channel<MapStatus> = ch_mapped
        .filter { r -> r.pass != null }
        .map { r -> record(id: r.id, pass: r.pass as Boolean) }

    def ch_percent_mapped_pass: Channel<PercentMappedPass> = ch_mapped.map { r -> record(id: r.id, percent_mapped: r.percent_mapped, pass: r.pass) }

    // Samples without a mapping percentage are never filtered
    ch_genome_bam = ch_mapped.filter { r -> r.pass == null || r.pass }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //

    // Some tools (Ex. Parabricks) may have already run marked duplicates during alignment
    def markdups_done = !params.skip_markduplicates && params.use_parabricks_star
    def ch_markdup: Channel<MarkdupBam> = channel.empty()
    if (!params.skip_markduplicates && !params.with_umi && !markdups_done) {
        ch_markdup = BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            ch_fasta,
            ch_fai,
            !params.use_rustqc
        )

        // Only bam, bai and metrics are merged: joining the whole result would overwrite the aligner's samtools stats with null when RustQC skips them
        ch_genome_marked = ch_genome_bam.join(
            ch_markdup.map { r -> record(id: r.id, bam: r.bam, bai: r.bai, metrics: r.metrics) },
            by: 'id',
            remainder: true
        )
        ch_genome_marked.subscribe { r ->
            if( r.meta == null || r.bam == null || r.metrics == null ) {
                error "Sample '${r.id}' is missing its BAM from mark duplicates"
            }
        }
        ch_genome_bam = ch_genome_marked.filter { r -> r.meta != null && r.bam != null && r.metrics != null }

        ch_mqc_files = ch_mqc_files.mix(
            ch_markdup
                .filter { r -> r.samtools != null }
                .map { r -> record(id: r.id, files: [r.samtools.stats, r.samtools.flagstat, r.samtools.idxstats, r.metrics]) }
        )
        ch_mqc_report_only = ch_mqc_report_only.mix(
            ch_markdup.filter { r -> r.samtools == null }.map { r -> r.metrics }
        )
    }

    //
    // MODULE: StringTie assembly and quantification
    //
    def ch_stringtie: Channel<StringtieSample>          = channel.empty()
    def ch_stringtie_merged: Channel<StringtieMerged>   = channel.empty()
    if (!params.skip_stringtie) {
        def ch_stringtie_input: Channel<StringtieInput> = ch_genome_bam.map { r -> r + record(lrbam: null) }

        if (params.stringtie_ignore_gtf) {
            stringtie_merge = BAM_STRINGTIE_MERGE(
                ch_stringtie_input,
                channel.value([]),
                ch_gtf
            )
            ch_stringtie_merged = stringtie_merge
            ch_stringtie_gtf = stringtie_merge.collect().map { merged -> merged.isEmpty() ? null : merged.toList()[0].merged_gtf }
        } else {
            ch_stringtie_gtf = ch_gtf
        }
        def ch_stringtie_samples: Channel<StringtieResult> = STRINGTIE_STRINGTIE(
            ch_stringtie_input,
            channel.value(['expression-estimation']),
            ch_stringtie_gtf
        )

        // Per-sample de novo assemblies that fed the merged GTF
        if (params.stringtie_ignore_gtf) {
            ch_stringtie = ch_stringtie_samples
                .join(stringtie_merge.flatMap { r -> r.assemblies }.map { r -> record(id: r.id, denovo: r) }, by: 'id')
        } else {
            ch_stringtie = ch_stringtie_samples.map { r -> r + record(denovo: null) }
        }
    }

    //
    // Pre-compute param-derived values for QC subworkflow
    //
    def biotype = (params.gencode ? "gene_type" : params.featurecounts_group_type) as String
    def rseqc_modules = qc_tools.findAll { tool -> tool.startsWith('rseqc_') }.collect { tool -> tool.replace('rseqc_', '') }

    def ch_bam_qc                                       = channel.empty()
    def ch_bam_qc_rustqc: Channel<RustqcResult>         = channel.empty()
    def ch_inferexperiment: Channel<InferExperimentLog> = channel.empty()

    if (!params.skip_qc) {
        if (params.use_rustqc) {
            //
            // MODULE: RustQC - single-pass replacement for multiple QC tools
            //
            ch_bam_qc_rustqc = RUSTQC (
                ch_genome_bam,
                ch_gtf,
            )

            // Drop non-MultiQC files. Excluding `*.featureCounts.tsv.summary`
            // keeps only the biotype summary, matching the default pipeline's
            // `featureCounts -g gene_biotype` output.
            ch_mqc_files = ch_mqc_files.mix(
                ch_bam_qc_rustqc.map { r ->
                    record(
                        id:    r.id,
                        files: r.all_files.findAll { f ->
                            !f.name.endsWith('.featureCounts.tsv.summary') &&
                                ((f.name =~ /(?i)\.(txt|tsv|xls|log|stats|flagstat|idxstats|html)$/).find() || f.name.contains('_mqc.'))
                        }.toList()
                    )
                }
            )

            // Extract infer_experiment from rseqc channel
            ch_inferexperiment = ch_bam_qc_rustqc
                .filter { r -> r.rseqc.inferexperiment != null }
                .map { r -> record(id: r.id, meta: r.meta, inferexperiment: r.rseqc.inferexperiment) }
        } else {
            //
            // SUBWORKFLOW: Post-alignment QC
            //
            ch_bam_qc = BAM_QC_RNASEQ (
                ch_genome_bam,
                ch_gtf,
                ch_gene_bed,
                ch_fasta,
                ch_fai,
                ch_biotypes_header_multiqc,
                qc_tools,
                biotype
            )

            ch_mqc_files = ch_mqc_files.mix(ch_bam_qc.map { r -> record(id: r.id, files: r.mqc_files) })
            ch_inferexperiment = ch_bam_qc
                .filter { r -> r.rseqc != null && r.rseqc.inferexperiment != null }
                .map { r -> record(id: r.id, meta: r.meta, inferexperiment: r.rseqc.inferexperiment) }
        }
    }

    //
    // Build the per-sample strand-classification record consumed by the
    // MultiQC Strandedness checks section. When RSeQC / RustQC ran we
    // classify via `classifyStrand`; otherwise we surface Salmon's
    // auto-inference so --skip_rseqc / --skip_qc users still see the
    // available signal. RustQC always emits infer_experiment regardless
    // of rseqc_modules.
    //
    def run_infer_experiment = !params.skip_qc && (params.use_rustqc || rseqc_modules.contains('infer_experiment'))
    def ch_strand_data: Channel<StrandData>     = channel.empty()
    def ch_strand_status: Channel<StrandStatus> = channel.empty()
    if (run_infer_experiment) {
        ch_strand_data = ch_inferexperiment.map { r ->
            classifyStrand(r.meta, r.inferexperiment, params.stranded_threshold, params.unstranded_threshold)
        }
        ch_strand_status = ch_strand_data.map { r -> record(id: r.id, pass: r.status == 'pass') }
    }
    else {
        ch_strand_data = ch_reads_ok
            .filter { r -> r.meta.salmon_strand_analysis != null }
            .map { r -> record(id: r.id, meta: r.meta, provided: 'auto', status: '-', salmon: r.meta.salmon_strand_analysis, rseqc: null) }
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    // Stranded libraries get per-strand + combined bigWigs; unstranded libraries get only the combined one.
    //
    def ch_bigwig: Channel<BigwigSample> = channel.empty()
    if (!params.skip_bigwig) {
        ch_genomecov_input = ch_genome_bam.map { r -> record(id: r.id, meta: r.meta, intervals: r.bam, scale: 1) }
        ch_genomecov_input_stranded = ch_genomecov_input.filter { r -> r.meta.strandedness in ['forward', 'reverse'] }

        ch_bedgraph_fw = BEDTOOLS_GENOMECOV_FW (
            ch_genomecov_input_stranded,
            null,
            'bedGraph',
            true
        )
        ch_bedgraph_rev = BEDTOOLS_GENOMECOV_REV (
            ch_genomecov_input_stranded,
            null,
            'bedGraph',
            true
        )
        ch_bedgraph_combined = BEDTOOLS_GENOMECOV_COMBINED (
            ch_genomecov_input,
            null,
            'bedGraph',
            true
        )

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        ch_bigwig_fw = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            ch_bedgraph_fw,
            ch_chrom_sizes
        )

        ch_bigwig_rev = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            ch_bedgraph_rev,
            ch_chrom_sizes
        )

        ch_bigwig_combined = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_COMBINED (
            ch_bedgraph_combined,
            ch_chrom_sizes
        )

        // Every sample gets a combined track; only stranded ones get forward/reverse
        ch_bigwig = ch_bigwig_combined
            .map { r -> record(id: r.id, meta: r.meta, combined: r, forward: null, reverse: null) }
            .join(ch_bigwig_fw.map { r -> record(id: r.id, forward: r) }, by: 'id', remainder: true)
            .join(ch_bigwig_rev.map { r -> record(id: r.id, reverse: r) }, by: 'id', remainder: true)
    }

    def ch_contaminants: Channel<Contaminants> = channel.empty()
    if (!params.skip_qc) {
        //
        // Contaminant screening (Kraken2/Bracken/Sylph)
        //
        if (params.contaminant_screening_input == 'trim_only') {
            ch_contaminant_reads = ch_preprocessed
                .filter { r -> r.reads_trimmed != null }
                .map { r -> record(id: r.id, meta: r.meta, reads: r.reads_trimmed) }
        } else if (params.contaminant_screening_input == 'trimmed') {
            ch_contaminant_reads = ch_reads_ok
        } else if (params.contaminant_screening_input == 'raw') {
            ch_contaminant_reads = ch_preprocessed.map { r -> record(id: r.id, meta: r.meta, reads: r.reads_cat) }
        } else {
            ch_contaminant_reads = ch_star.mix(ch_hisat2)
                .filter { r -> !r.unmapped.isEmpty() }
                .map { r -> record(id: r.id, meta: r.meta, reads: r.unmapped) }
        }

        if (params.contaminant_screening in ['kraken2', 'kraken2_bracken'] ) {
            ch_kraken2 = KRAKEN2 (
                ch_contaminant_reads,
                ch_kraken_db,
                params.save_kraken_assignments,
                params.save_kraken_unassigned
            )

            ch_contaminants = ch_kraken2.map { r ->
                record(id: r.id, meta: r.meta, kraken2: r, bracken: null, sylph: null, sylphtax: null)
            }

            if (params.contaminant_screening == 'kraken2') {
                ch_mqc_files = ch_mqc_files.mix(ch_kraken2.map { r -> record(id: r.id, files: [r.report]) })
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                ch_bracken = BRACKEN (
                    ch_kraken2,
                    ch_kraken_db
                )
                ch_mqc_files = ch_mqc_files.mix(ch_bracken.map { r -> record(id: r.id, files: [r.report]) })
                ch_contaminants = ch_contaminants
                    .join(ch_bracken.map { r -> record(id: r.id, bracken: r) }, by: 'id')
            }
        } else if (params.contaminant_screening == 'sylph') {
            def sylph_databases = (params.sylph_db ? params.sylph_db.split(',').collect{ path -> file(path.trim()) } : []) as List<Path>
            ch_sylph_databases = channel.value(sylph_databases)
            ch_sylph = SYLPH_PROFILE (
                ch_contaminant_reads,
                ch_sylph_databases
            )
            ch_sylph_profile = ch_sylph.filter{ r -> !r.profile_out.isEmpty() }

            def sylph_taxonomies = (params.sylph_taxonomy ? params.sylph_taxonomy.split(',').collect{ path -> file(path.trim()) } : []) as List<Path>
            ch_sylph_taxonomies = channel.value(sylph_taxonomies)
            ch_sylphtax = SYLPHTAX_TAXPROF (
                ch_sylph_profile,
                ch_sylph_taxonomies
            )
            ch_mqc_files = ch_mqc_files.mix(ch_sylphtax.map { r -> record(id: r.id, files: [r.taxprof_output]) })

            // Every profile is published, but empty ones never reach SYLPHTAX_TAXPROF
            ch_contaminants = ch_sylph
                .map { r -> record(id: r.id, meta: r.meta, kraken2: null, bracken: null, sylph: record(profile: r.profile_out), sylphtax: null) }
                .join(ch_sylphtax.map { r -> record(id: r.id, sylphtax: record(taxprof: r.taxprof_output)) }, by: 'id', remainder: true)
        }
    }

    //
    // SUBWORKFLOW: Pseudoalignment and quantification with Salmon
    //
    def ch_quant_pseudo: Channel<SalmonQuantSample>        = channel.empty()
    def ch_quant_pseudo_kallisto: Channel<KallistoQuantSample> = channel.empty()
    def ch_quant_merged_pseudo: Channel<QuantMerged>       = channel.empty()
    def ch_deseq2_pseudo: Channel<Deseq2Qc>                = channel.empty()
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {

        if (params.pseudo_aligner == 'salmon') {
            ch_pseudo_index = ch_salmon_index
        } else {
            ch_pseudo_index = ch_kallisto_index
        }

        pseudo = QUANTIFY_PSEUDO_ALIGNMENT (
            ch_samplesheet,
            ch_reads_ok,
            ch_pseudo_index,
            null,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.pseudo_aligner,
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge
        )

        ch_quant_pseudo          = pseudo.salmon
        ch_quant_pseudo_kallisto = pseudo.kallisto
        ch_quant_merged_pseudo   = pseudo.merged

        // MultiQC parses the Salmon quant directory and the Kallisto log
        def ch_pseudo_mqc: Channel<MultiqcFiles> = ch_quant_pseudo
            .map { r -> record(id: r.id, files: [r.quant_dir]) }
            .mix(ch_quant_pseudo_kallisto.map { r -> record(id: r.id, files: [r.log]) })
        ch_mqc_files = ch_mqc_files.mix(ch_pseudo_mqc)

        if (run_deseq2_qc) {
            ch_deseq2_pseudo = DESEQ2_QC_PSEUDO (
                ch_quant_merged_pseudo,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
        }
    }

    ch_mqc_report_only = ch_mqc_report_only.mix(
        ch_deseq2.mix(ch_deseq2_pseudo).flatMap { r -> [r.pca_multiqc, r.dists_multiqc].findAll { f -> f != null } }
    )

    //
    // Collate and save software versions from the `versions` topic. Entries are either
    // `path(versions.yml)` (legacy file-emit style) or `(task.process, tool, version)`
    // tuples (inline `eval` style).
    //
    ch_topic_versions = channel.topic('versions').unique()

    ch_versions_string = ch_topic_versions
        .filter { entry -> !(entry instanceof Path) }
        .map { entry -> tuple(entry[0].tokenize(':').last(), "  ${entry[1]}: ${entry[2]}") }
        .groupBy()
        .map { process, tool_versions ->
            "${process}:\n${tool_versions.toSet().toSorted().join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_topic_versions.filter { entry -> entry instanceof Path })
        .mix(ch_versions_string)
        .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .map { p -> p as Path }

    def ch_pipeline_info: Channel<PipelineInfo> = ch_collated_versions.map { versions -> record(versions: versions) }

    //
    // SUBWORKFLOW: MultiQC
    //
    def ch_multiqc: Channel<MultiqcReport> = channel.empty()
    def ch_multiqc_report: Channel<Path>   = channel.empty()

    if (!params.skip_multiqc) {
        ch_multiqc = MULTIQC_RNASEQ(
            ch_input.map { s -> s.id },
            ch_mqc_files,
            ch_mqc_sample_only,
            ch_mqc_report_only,
            ch_strand_data,
            ch_trim_read_count,
            ch_percent_mapped_pass,
            aligner_display_name,
            ch_fastq,
            ch_collated_versions,
            file("$projectDir/assets/multiqc_config.yml", checkIfExists: true),
            params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : null,
            params.multiqc_logo   ? file(params.multiqc_logo,   checkIfExists: true) : null,
            params.multiqc_methods_description
                ? file(params.multiqc_methods_description)
                : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true),
            file("$projectDir/assets/strand_check_summary.yaml",     checkIfExists: true),
            file("$projectDir/assets/strand_check_composition.yaml", checkIfExists: true),
            sample_status_header_multiqc,
            params.min_trimmed_reads,
            params.skip_quantification_merge
        )
        ch_multiqc_report = ch_multiqc.map { r -> r.report }
    }

    emit:
    trim_status:           Channel<TrimStatus>                       = ch_trim_status           // pass: meets min_trimmed_reads
    map_status:            Channel<MapStatus>                        = ch_map_status            // pass: meets min_mapped_reads; samples with a mapping percentage only
    strand_status:         Channel<StrandStatus>                     = ch_strand_status         // pass: strandedness check passed
    multiqc_report:        Channel<Path>                             = ch_multiqc_report
    reads:                 Channel<SampleRuns>                       = ch_fastq
    percent_mapped:        Channel<PercentMapped>                    = ch_percent_mapped

    // Stage result records, keyed on id
    preprocessed:          Channel<FastqQcTrimFilterSetstrandedness> = ch_preprocessed
    aligned:               Channel<AlignedSample>                    = ch_aligned               // StarAligned | Bowtie2Aligned | Hisat2Aligned
    umi_dedup:             Channel<UmiDedupBam>                      = ch_umi_dedup
    markdup:               Channel<MarkdupBam>                       = ch_markdup
    bam_qc:                Channel<BamQcRnaseq>                      = ch_bam_qc
    bam_qc_rustqc:         Channel<RustqcResult>                     = ch_bam_qc_rustqc
    quant:                 Channel<RsemQuantSample>                  = ch_quant
    quant_salmon:          Channel<SalmonQuantSample>                = ch_quant_salmon          // Salmon on the transcriptome BAM
    quant_merged:          Channel<QuantMerged>                      = ch_quant_merged          // alignment-based quantifier
    quant_rsem_merge:      Channel<RsemMergeSample>                  = ch_quant_rsem_merge      // empty unless --aligner star_rsem
    quant_pseudo:          Channel<SalmonQuantSample>                = ch_quant_pseudo          // Salmon pseudo-alignment
    quant_pseudo_kallisto: Channel<KallistoQuantSample>              = ch_quant_pseudo_kallisto // Kallisto pseudo-alignment
    quant_merged_pseudo:   Channel<QuantMerged>                      = ch_quant_merged_pseudo   // pseudo-aligner
    contaminants:          Channel<Contaminants>                     = ch_contaminants
    stringtie:             Channel<StringtieSample>                  = ch_stringtie
    bigwig:                Channel<BigwigSample>                     = ch_bigwig

    // Run-level result records
    stringtie_merged:      Channel<StringtieMerged>                  = ch_stringtie_merged
    deseq2:                Channel<Deseq2Qc>                         = ch_deseq2                // alignment-based quantifier
    deseq2_pseudo:         Channel<Deseq2Qc>                         = ch_deseq2_pseudo         // pseudo-aligner
    rrna_references:       Value<RrnaReferences>                     = ch_rrna_references
    multiqc:               Channel<MultiqcReport>                    = ch_multiqc               // per sample under skip_quantification_merge
    pipeline_info:         Channel<PipelineInfo>                     = ch_pipeline_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
