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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
params.bowtie2_index    = getGenomeAttribute('bowtie2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                     } from './workflows/rnaseq'
include { PREPARE_GENOME_REFERENCES  } from './subworkflows/local/prepare_genome_references'
include { PREPARE_GENOME_INDICES     } from './subworkflows/local/prepare_genome_indices'
include { PIPELINE_INITIALISATION    } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION        } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { defineQcTools              } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy          } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

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

    //
    // SUBWORKFLOW: Prepare reference genome files (FASTA, GTF, BED, transcript FASTA, chrom.sizes, rRNA FASTAs, Kraken DB)
    //
    PREPARE_GENOME_REFERENCES (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.ribo_database_manifest,
        params.kraken_db,
        params.gencode,
        params.gffread_transcript_fasta,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        params.remove_ribo_rna && !(params.ribo_removal_tool == "bowtie2" && params.bowtie2_rrna_index) ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.contaminant_screening,
        params.prokaryotic ?: false
    )

    //
    // SUBWORKFLOW: Build or load aligner / pseudo-aligner / filtering indices
    //
    PREPARE_GENOME_INDICES (
        PREPARE_GENOME_REFERENCES.out.fasta_fai,
        PREPARE_GENOME_REFERENCES.out.gtf,
        PREPARE_GENOME_REFERENCES.out.transcript_fasta,
        PREPARE_GENOME_REFERENCES.out.rrna_fastas,
        params.fasta ? true : false,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bowtie2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.bowtie2_rrna_index,
        params.aligner,
        params.pseudo_aligner,
        params.skip_bbsplit,
        params.remove_ribo_rna ? params.ribo_removal_tool : null,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.use_parabricks_star,
        isStarIndexLegacy() ?: false,
        params.hisat2_build_memory
    )

    ch_genome = PREPARE_GENOME_REFERENCES.out.results
        .combine(PREPARE_GENOME_INDICES.out.results)
        .map { references, indices -> references + record(index: indices) }

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME_REFERENCES
            .out
            .fasta_fai
            .map { _meta, _fasta, fai -> checkMaxContigSize(fai) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    def qc_tools = defineQcTools(params)

    RNASEQ (
        ch_samplesheet,
        PREPARE_GENOME_REFERENCES.out.fasta_fai,
        PREPARE_GENOME_REFERENCES.out.gtf,
        PREPARE_GENOME_REFERENCES.out.chrom_sizes,
        PREPARE_GENOME_REFERENCES.out.gene_bed,
        PREPARE_GENOME_REFERENCES.out.transcript_fasta,
        PREPARE_GENOME_INDICES.out.star_index,
        PREPARE_GENOME_INDICES.out.rsem_index,
        PREPARE_GENOME_INDICES.out.hisat2_index,
        PREPARE_GENOME_INDICES.out.bowtie2_index,
        PREPARE_GENOME_INDICES.out.salmon_index,
        PREPARE_GENOME_INDICES.out.kallisto_index,
        PREPARE_GENOME_INDICES.out.bbsplit_index,
        PREPARE_GENOME_REFERENCES.out.rrna_fastas,
        PREPARE_GENOME_INDICES.out.sortmerna_index,
        PREPARE_GENOME_INDICES.out.bowtie2_rrna_index,
        PREPARE_GENOME_INDICES.out.splicesites,
        PREPARE_GENOME_REFERENCES.out.kraken_db,
        qc_tools
    )

    ch_genome = ch_genome
        .combine(RNASEQ.out.rrna_references)
        .map { genome, rrna_references -> genome + record(rrna_references: rrna_references) }

    // Indexes built inside preprocessing share their directory name with the genome-level ones
    // (`idx`, `salmon`), so they cannot ride the genome target (nextflow-io/nextflow#6617).
    ch_genome_preprocessing = RNASEQ.out.preprocessing_references
        .filter { r -> r.salmon_index != null || r.sortmerna_index != null }

    // Same-basename fields split out per stage to avoid a >> rename-key collision (nextflow-io/nextflow#6617).
    ch_lint_raw     = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.raw) }.filter { s -> s.file != null }
    ch_lint_trimmed = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.trimmed) }.filter { s -> s.file != null }
    ch_lint_bbsplit = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.bbsplit) }.filter { s -> s.file != null }
    ch_lint_ribo    = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.ribo) }.filter { s -> s.file != null }

    // CUSTOM_RSEMMERGECOUNTS and tximport both write rsem.merged.* basenames. A target publishes every
    // file in its records and >> keys the rename on basename (nextflow-io/nextflow#6617), so the
    // quant_merged records must not carry rsem_merge at all; it gets its own target.
    ch_quant_rsem_merge = RNASEQ.out.quant_merged
        .filter { r -> r.rsem_merge != null }
        .map { r -> record(id: r.id, rsem_merge: r.rsem_merge) }
    ch_quant_merged = RNASEQ.out.quant_merged.map { r -> r + record(rsem_merge: null) }

    // The prepared rRNA FASTAs publish whether or not --save_reference is set, so they cannot ride the genome target.
    ch_rrna_seqkit = RNASEQ.out.rrna_references.filter { r -> r.seqkit_prefixed || r.seqkit_converted }

    // Flattened to one file per record so each can carry its own precomputed rename-form >> target.
    ch_bam_qc_rustqc_files = RNASEQ.out.bam_qc_rustqc
        .flatMap { s ->
            (s.samtools ?: []).collect      { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'samtools', f)) } +
            (s.dupradar ?: []).collect      { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'dupradar', f)) } +
            (s.featurecounts ?: []).collect { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'featurecounts', f)) } +
            (s.preseq ?: []).collect        { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'preseq', f)) } +
            (s.rseqc ?: []).collect         { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'rseqc', f)) } +
            (s.qualimap ?: []).collect      { f -> record(id: s.id, file: f, target: rustqcTarget(s, 'qualimap', f)) }
        }
        .filter { r -> r.target != null }

    // samplesheet_with_bams.csv rows: one per sequencing run of a sample, with the aligned record's
    // meta so an inferred strandedness replaces 'auto'. genome_bam is the
    // coordinate-sorted BAM; bowtie2_salmon aligns to the transcriptome, so its unsorted bowtie2 BAM is the transcriptome_bam.
    // Filtered to samples that actually went through alignment here (excludes
    // the BAM-input passthrough placeholder record, which carries no orig_bam).
    ch_samplesheet_rows = RNASEQ.out.aligned
        .filter { r -> r.orig_bam != null }
        .map { r -> [r.id, r] }
        .join(RNASEQ.out.reads.map { meta, runs -> [meta.id, meta, runs] })
        .join(RNASEQ.out.percent_mapped)
        .flatMap { sample_id, r, _meta, runs, percent_mapped ->
            runs.collect { run ->
                record(
                    sample:            sample_id,
                    fastq_1:           run[0],
                    fastq_2:           run.size() > 1 ? run[1] : null,
                    strandedness:      r.meta.strandedness,
                    seq_platform:      r.meta.seq_platform ?: params.seq_platform,
                    seq_center:        r.meta.seq_center ?: params.seq_center,
                    genome_bam:        r.bam,
                    percent_mapped:    percent_mapped,
                    transcriptome_bam: params.aligner == 'bowtie2_salmon' ? r.orig_bam : r.transcriptome_bam
                )
            }
        }

    emit:
    trim_status         = RNASEQ.out.trim_status         // channel: [id, boolean]
    map_status          = RNASEQ.out.map_status          // channel: [id, boolean]
    strand_status       = RNASEQ.out.strand_status       // channel: [id, boolean]
    multiqc_report      = RNASEQ.out.multiqc_report      // channel: /path/to/multiqc_report.html
    genome              = ch_genome                      // channel: GenomeReferences fields + index: GenomeIndices + rrna_references
    genome_preprocessing = ch_genome_preprocessing       // channel: record(salmon_index, sortmerna_index), only when preprocessing built at least one
    rrna_seqkit         = ch_rrna_seqkit                 // channel: record(bowtie2_index, seqkit_prefixed, seqkit_converted), only when the bowtie2 rRNA index is built

    // Stage result records, keyed on id
    preprocessed        = RNASEQ.out.preprocessed        // channel: FastqQcTrimFilterSetstrandedness
    lint_raw            = ch_lint_raw                    // channel: record(id, file), FQ_LINT on raw reads
    lint_trimmed        = ch_lint_trimmed                // channel: record(id, file), FQ_LINT on trimmed reads
    lint_bbsplit        = ch_lint_bbsplit                // channel: record(id, file), FQ_LINT on BBSplit-filtered reads
    lint_ribo           = ch_lint_ribo                   // channel: record(id, file), FQ_LINT on rRNA-removed reads
    aligned             = RNASEQ.out.aligned             // channel: StarAligned | Bowtie2Aligned | Hisat2Aligned
    umi_dedup           = RNASEQ.out.umi_dedup           // channel: UmiDedupBam
    markdup             = RNASEQ.out.markdup             // channel: MarkdupBam
    bam_qc              = RNASEQ.out.bam_qc              // channel: BamQcRnaseq
    bam_qc_rustqc       = ch_bam_qc_rustqc_files         // channel: record(id, file, target), one entry per RustQC output file
    samplesheet         = ch_samplesheet_rows            // channel: record(sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam), one entry per sequencing run
    quant               = RNASEQ.out.quant               // channel: RsemQuantSample | PseudoQuantSample, alignment-based quantifier
    quant_merged        = ch_quant_merged                // channel: RsemQuantMerged | QuantMerged with rsem_merge null, alignment-based quantifier
    quant_rsem_merge    = ch_quant_rsem_merge            // channel: record(id, rsem_merge: RsemMerge), CUSTOM_RSEMMERGECOUNTS outputs
    quant_pseudo        = RNASEQ.out.quant_pseudo        // channel: PseudoQuantSample, pseudo-aligner
    quant_merged_pseudo = RNASEQ.out.quant_merged_pseudo // channel: QuantMerged, pseudo-aligner
    contaminants        = RNASEQ.out.contaminants        // channel: record(id, meta, kraken2, bracken, sylph, sylphtax)
    stringtie           = RNASEQ.out.stringtie           // channel: record(id, meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo: StringtieAssembly?)
    bigwig              = RNASEQ.out.bigwig              // channel: record(id, meta, combined, forward, reverse), each BigwigFiles

    // Run-level result records
    stringtie_merged    = RNASEQ.out.stringtie_merged    // channel: StringtieMerged
    deseq2              = RNASEQ.out.deseq2              // channel: record(rdata, pca_vals, plots_pdf, sample_dists, size_factors, log), alignment-based quantifier
    deseq2_pseudo       = RNASEQ.out.deseq2_pseudo       // channel: record(rdata, pca_vals, plots_pdf, sample_dists, size_factors, log), pseudo-aligner
    multiqc             = RNASEQ.out.multiqc             // channel: MultiqcReport, per sample under skip_quantification_merge
    pipeline_info       = RNASEQ.out.pipeline_info       // channel: record(versions)
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
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
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
        NFCORE_RNASEQ.out.multiqc_report,
        NFCORE_RNASEQ.out.trim_status,
        NFCORE_RNASEQ.out.map_status,
        NFCORE_RNASEQ.out.strand_status
    )

    publish:
    contaminants     = NFCORE_RNASEQ.out.contaminants
    stringtie        = NFCORE_RNASEQ.out.stringtie
    stringtie_merged = NFCORE_RNASEQ.out.stringtie_merged
    bigwig           = NFCORE_RNASEQ.out.bigwig
    genome           = NFCORE_RNASEQ.out.genome
    genome_preprocessing = NFCORE_RNASEQ.out.genome_preprocessing
    rrna_seqkit      = NFCORE_RNASEQ.out.rrna_seqkit
    preprocessed     = NFCORE_RNASEQ.out.preprocessed
    lint_raw         = NFCORE_RNASEQ.out.lint_raw
    lint_trimmed     = NFCORE_RNASEQ.out.lint_trimmed
    lint_bbsplit     = NFCORE_RNASEQ.out.lint_bbsplit
    lint_ribo        = NFCORE_RNASEQ.out.lint_ribo
    aligned          = NFCORE_RNASEQ.out.aligned
    umi_dedup        = NFCORE_RNASEQ.out.umi_dedup
    markdup          = NFCORE_RNASEQ.out.markdup
    samplesheet      = NFCORE_RNASEQ.out.samplesheet
    quant               = NFCORE_RNASEQ.out.quant
    quant_merged        = NFCORE_RNASEQ.out.quant_merged
    quant_rsem_merge    = NFCORE_RNASEQ.out.quant_rsem_merge
    quant_pseudo        = NFCORE_RNASEQ.out.quant_pseudo
    quant_merged_pseudo = NFCORE_RNASEQ.out.quant_merged_pseudo
    deseq2              = NFCORE_RNASEQ.out.deseq2
    deseq2_pseudo       = NFCORE_RNASEQ.out.deseq2_pseudo
    bam_qc              = NFCORE_RNASEQ.out.bam_qc
    bam_qc_rustqc       = NFCORE_RNASEQ.out.bam_qc_rustqc
    multiqc             = NFCORE_RNASEQ.out.multiqc
    pipeline_info       = NFCORE_RNASEQ.out.pipeline_info
}

// Run-level records (e.g. a cross-sample merged file) are never sample-prefixed.
def samplePrefix(r) { params.skip_quantification_merge ? "${r.id}/" : '' }
def alignerDir(r)    { "${samplePrefix(r)}${params.aligner}" }
def trimLogDir(s)    { params.trimmer == 'fastp' ? "${samplePrefix(s)}${params.trimmer}/log" : "${samplePrefix(s)}${params.trimmer}" }
def saveAlignBam(_s)  { params.save_align_intermeds || params.skip_markduplicates }
def saveUmiBam(_s)    { params.save_align_intermeds || params.save_umi_intermeds }
def umiDedupToolDir(_s) { params.umi_dedup_tool == 'umicollapse' ? 'umicollapse' : 'umitools' }
def pseudoAlignerDir(r) { "${samplePrefix(r)}${params.pseudo_aligner}" }

def multiqcDir(m) {
    def suffix = params.skip_alignment ? '' : "/${params.aligner}"
    m.id == 'multiqc_report' ? "multiqc${suffix}" : "${m.id}/multiqc${suffix}"
}

// True whenever at least one field of a `preprocessed` record is guaranteed
// non-null: the disjunction of the conditions on the `preprocessed` path
// closure's `>>` lines. Mirrors the skip/enable args passed into
// FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS in workflows/rnaseq/main.nf.
def preprocessedPublishes() {
    !(params.skip_fastqc || params.skip_qc) ||
    !params.skip_trimming ||
    params.with_umi ||
    (!params.skip_bbsplit && params.fasta) ||
    params.remove_ribo_rna ||
    params.save_merged_fastq
}

// Directory for preprocessed.reads when rRNA removal produced them, or null.
// Trimmed, UMI-extracted and BBSplit reads publish from their own fields.
def rrnaFilteredReadsDir(s) {
    if (!params.save_non_ribo_reads)      { return null }
    if (s.rrna?.sortmerna_log != null)    { return "${samplePrefix(s)}sortmerna/" }
    if (s.rrna?.ribodetector_log != null) { return "${samplePrefix(s)}ribodetector/" }
    // Only the single-end --un-gz FASTQs are published; paired-end reads rebuilt by SAMTOOLS_FASTQ_BOWTIE2 are not.
    if (s.rrna?.bowtie2_log != null)      { return s.meta.single_end ? "${samplePrefix(s)}bowtie2_rrna/" : null }
    return null
}

// Path of a RustQC output file below its <id>/<category>/ task directory.
def rustqcRelPath(id, category, file) {
    def marker = "${id}/${category}/"
    def path = file.toString()
    path.substring(path.lastIndexOf(marker) + marker.length())
}

// Computes the per-tool destination directory for one RustQC output file.
def rustqcTarget(r, category, file) {
    def dir = "${alignerDir(r)}/rustqc"
    def base = file.name
    if (category == 'samtools')     { return "${dir}/samtools_stats/${base}" }
    if (category == 'preseq')       { return "${dir}/preseq/${base}" }
    if (category == 'dupradar') {
        if (base.contains('Boxplot'))                                { return "${dir}/dupradar/box_plot/${base}" }
        if (base.contains('ExpDens') && !base.contains('Curve_mqc')) { return "${dir}/dupradar/scatter_plot/${base}" }
        if (base.contains('expressionHist'))                        { return "${dir}/dupradar/histogram/${base}" }
        if (base.contains('dupMatrix'))                              { return "${dir}/dupradar/gene_data/${base}" }
        if (base.contains('intercept_slope'))                        { return "${dir}/dupradar/intercepts_slope/${base}" }
        return "${dir}/dupradar/${base}"
    }
    if (category == 'featurecounts') {
        if (base.endsWith('.featureCounts.biotype.tsv.summary')) { return "${dir}/featurecounts/${base.replace('.biotype.tsv.summary', '.tsv.summary')}" }
        if (base.endsWith('.featureCounts.tsv.summary'))         { return null }
        return "${dir}/featurecounts/${base}"
    }
    if (category == 'rseqc') {
        def relPath = rustqcRelPath(r.id, category, file)
        def tool = relPath.tokenize('/')[0]
        if (tool in ['junction_annotation', 'junction_saturation', 'inner_distance', 'read_duplication']) {
            if (base.endsWith('.r'))                                { return "${dir}/rseqc/${tool}/rscript/${base}" }
            if (base.endsWith('.png') || base.endsWith('.svg'))     { return "${dir}/rseqc/${tool}/plot/${base}" }
            if (base.endsWith('.xls'))                              { return "${dir}/rseqc/${tool}/xls/${base}" }
            if (base.endsWith('.bed'))                              { return "${dir}/rseqc/${tool}/bed/${base}" }
            if (base.endsWith('.junction_annotation.log'))          { return "${dir}/rseqc/${tool}/log/${base}" }
            if (base.endsWith('.txt'))                              { return "${dir}/rseqc/${tool}/txt/${base}" }
        }
        return "${dir}/rseqc/${relPath}"
    }
    if (category == 'qualimap') {
        return "${dir}/qualimap/${r.id}/${rustqcRelPath(r.id, category, file)}"
    }
    return null
}

output {
    contaminants {   // record(id, meta, kraken2, bracken, sylph, sylphtax); exactly one tool branch is populated per run
        enabled !params.skip_qc && params.contaminant_screening
        path { s ->
            s.kraken2?.report                      >> "${alignerDir(s)}/contaminants/kraken2/kraken_reports/"
            s.kraken2?.classified_reads_fastq      >> (params.save_kraken_assignments ? "${alignerDir(s)}/contaminants/kraken2/kraken_reports/" : null)
            s.kraken2?.unclassified_reads_fastq    >> (params.save_kraken_assignments ? "${alignerDir(s)}/contaminants/kraken2/kraken_reports/" : null)
            s.kraken2?.classified_reads_assignment >> (params.save_kraken_unassigned ? "${alignerDir(s)}/contaminants/kraken2/kraken_reports/" : null)
            s.bracken?.abundance                   >> "${alignerDir(s)}/contaminants/bracken/"
            s.bracken?.report                      >> "${alignerDir(s)}/contaminants/bracken/"
            s.sylph?.profile                       >> "${alignerDir(s)}/contaminants/sylph/"
            s.sylphtax?.taxprof                    >> "${alignerDir(s)}/contaminants/sylph/"
        }
    }

    stringtie {   // record(id, meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo: StringtieAssembly?)
        enabled !params.skip_stringtie
        path { s ->
            s.transcript_gtf >> "${alignerDir(s)}/stringtie/"
            s.abundance >> "${alignerDir(s)}/stringtie/"
            s.coverage_gtf >> "${alignerDir(s)}/stringtie/"
            s.ballgown >> "${alignerDir(s)}/stringtie/"
            s.denovo?.transcript_gtf >> "${alignerDir(s)}/stringtie/"
        }
    }

    stringtie_merged {   // record(id, merged_gtf); cross-sample, only under --stringtie_ignore_gtf, never sample-prefixed
        enabled !params.skip_stringtie && params.stringtie_ignore_gtf
        path { s -> s.merged_gtf >> "${params.aligner}/stringtie/" }
    }

    bigwig {   // record(id, meta, combined, forward, reverse), each BigwigFiles { bigwig, bedgraph }; bedgraph stays unrouted
        enabled !params.skip_bigwig
        path { s ->
            s.combined?.bigwig >> "${alignerDir(s)}/bigwig/"
            s.forward?.bigwig  >> "${alignerDir(s)}/bigwig/"
            s.reverse?.bigwig  >> "${alignerDir(s)}/bigwig/"
        }
    }

    genome {   // GenomeReferences + index: GenomeIndices + rrna_references; never sample-prefixed
        enabled params.save_reference
        path { g ->
            g.fasta >> 'genome/'
            g.fai >> 'genome/'
            g.gtf >> 'genome/'
            g.gene_bed >> 'genome/'
            g.transcript_fasta >> 'genome/'
            g.chrom_sizes >> 'genome/'
            g.rrna_fastas >> 'genome/'
            g.kraken_db >> 'genome/index/'
            g.intermediates?.gff >> 'genome/'
            g.intermediates?.additional_fasta >> 'genome/'
            g.intermediates?.gtf_pre_filter >> 'genome/'
            g.intermediates?.fasta_pre_concat >> 'genome/'
            g.intermediates?.gtf_pre_concat >> 'genome/'
            g.intermediates?.transcript_fasta_pre_gencode >> 'genome/'
            g.intermediates?.transcript_fasta_rsem_dir >> 'genome/'
            g.index?.star >> 'genome/index/'
            g.index?.rsem >> 'genome/index/'
            g.index?.rsem_transcript_fasta >> 'genome/index/'
            g.index?.hisat2 >> 'genome/index/'
            g.index?.hisat2_splicesites >> 'genome/index/'
            g.index?.bowtie2 >> 'genome/index/'
            g.index?.salmon >> 'genome/index/'
            g.index?.kallisto >> 'genome/index/'
            g.index?.bbsplit >> 'genome/index/'
            g.index?.bbsplit_log >> 'genome/index/'
            g.index?.sortmerna >> 'genome/sortmerna/'
            g.index?.bowtie2_rrna >> 'genome/index/'
            g.rrna_references?.bowtie2_index >> 'bowtie2_rrna/index/'
        }
    }

    genome_preprocessing {   // record(salmon_index, sortmerna_index); ch_genome_preprocessing guarantees at least one is non-null
        enabled params.save_reference
        path { p ->
            p.salmon_index >> 'genome/index/'
            p.sortmerna_index >> 'genome/sortmerna/'
        }
    }

    rrna_seqkit {   // record(bowtie2_index, seqkit_prefixed, seqkit_converted); ch_rrna_seqkit guarantees at least one FASTA list is non-empty
        path { r ->
            r.seqkit_prefixed >> 'seqkit/'
            r.seqkit_converted >> 'seqkit/'
        }
    }

    preprocessed {   // FastqQcTrimFilterSetstrandedness; no single anchor field survives every skip combination
        enabled preprocessedPublishes()
        path { s ->
            s.fastqc?.raw_html >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.raw_zip >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.trim_html >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.trim_zip >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.filtered_html >> "${samplePrefix(s)}fastqc/filtered/"
            s.fastqc?.filtered_zip >> "${samplePrefix(s)}fastqc/filtered/"
            s.trim?.html >> "${samplePrefix(s)}${params.trimmer}/"
            s.trim?.log >> "${trimLogDir(s)}/"
            s.trim?.json >> (params.trimmer == 'fastp' ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.trim?.unpaired >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.trim?.reads_fail >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.trim?.reads_merged >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.reads_trimmed >> (!params.skip_trimming && params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.umi?.log >> "${samplePrefix(s)}umitools/"
            s.umi?.reads >> (params.save_umi_intermeds ? "${samplePrefix(s)}umitools/" : null)
            s.bbsplit?.stats >> "${samplePrefix(s)}bbsplit/"
            s.bbsplit?.primary_reads >> (params.save_bbsplit_reads ? "${samplePrefix(s)}bbsplit/" : null)
            s.bbsplit?.other_genome_reads >> (params.save_bbsplit_reads ? "${samplePrefix(s)}bbsplit/" : null)
            s.rrna?.sortmerna_log >> "${samplePrefix(s)}sortmerna/"
            s.rrna?.ribodetector_log >> "${samplePrefix(s)}ribodetector/"
            s.rrna?.seqkit_stats >> "${samplePrefix(s)}ribodetector/"
            s.rrna?.bowtie2_log >> "${samplePrefix(s)}bowtie2_rrna/"
            s.reads_cat >> (params.save_merged_fastq ? "${samplePrefix(s)}fastq/" : null)
            s.reads >> rrnaFilteredReadsDir(s)
        }
    }

    // One target per FQ_LINT stage: same basename, avoids the >> collision (nextflow-io/nextflow#6617).
    lint_raw {   // record(id, file); file is guaranteed non-null, ch_lint_raw filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s)}fq_lint/raw/" }
    }

    lint_trimmed {   // record(id, file); file is guaranteed non-null, ch_lint_trimmed filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s)}fq_lint/trimmed/" }
    }

    lint_bbsplit {   // record(id, file); file is guaranteed non-null, ch_lint_bbsplit filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s)}fq_lint/bbsplit/" }
    }

    lint_ribo {   // record(id, file); file is guaranteed non-null, ch_lint_ribo filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s)}fq_lint/${params.ribo_removal_tool ?: 'sortmerna'}/" }
    }

    aligned {   // StarAligned | Bowtie2Aligned | Hisat2Aligned; anchor: samtools.stats
        path { s ->
            s.samtools?.stats >> "${alignerDir(s)}/samtools_stats/"
            s.samtools?.flagstat >> "${alignerDir(s)}/samtools_stats/"
            s.samtools?.idxstats >> "${alignerDir(s)}/samtools_stats/"
            s.bam >> (saveAlignBam(s) ? "${alignerDir(s)}/" : null)
            s.bai >> (saveAlignBam(s) ? "${alignerDir(s)}/" : null)
            s.orig_bam >> (params.save_align_intermeds ? "${alignerDir(s)}/" : null)
            s.transcriptome_bam >> (params.save_align_intermeds ? "${alignerDir(s)}/" : null)
            s.unmapped >> (params.save_unaligned ? "${alignerDir(s)}/unmapped/" : null)
            s.star?.log_final >> "${alignerDir(s)}/log/"
            s.star?.log_out >> "${alignerDir(s)}/log/"
            s.star?.log_progress >> "${alignerDir(s)}/log/"
            s.star?.tab >> "${alignerDir(s)}/log/"
            s.hisat2?.summary >> "${alignerDir(s)}/log/"
            s.bowtie2?.log >> "${alignerDir(s)}/log/"
            s.preexisting_bai >> "${samplePrefix(s)}samtools/"
        }
    }

    samplesheet {   // record(sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam); field order is the CSV column order
        enabled params.save_align_intermeds && !params.skip_alignment
        path { r ->
            r.genome_bam >> "${params.skip_quantification_merge ? "${r.sample}/" : ''}${params.aligner}/"
            r.transcriptome_bam >> "${params.skip_quantification_merge ? "${r.sample}/" : ''}${params.aligner}/"
        }
        index {
            path 'samplesheets/samplesheet_with_bams.csv'
            header true
        }
    }

    umi_dedup {   // UmiDedupBam; anchor: genome.stats
        path { s ->
            s.genome?.stats >> "${alignerDir(s)}/samtools_stats/"
            s.genome?.flagstat >> "${alignerDir(s)}/samtools_stats/"
            s.genome?.idxstats >> "${alignerDir(s)}/samtools_stats/"
            s.bam >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.bai >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.genomic_dedup_log >> "${alignerDir(s)}/${umiDedupToolDir(s)}/genomic_dedup_log/"
            s.tsv?.edit_distance >> "${alignerDir(s)}/umitools/"
            s.tsv?.per_umi >> "${alignerDir(s)}/umitools/"
            s.tsv?.umi_per_position >> "${alignerDir(s)}/umitools/"
            s.transcriptome?.coord_sorted_bam >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.transcriptome?.coord_sorted_bam_index >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.transcriptome?.coord_sorted_samtools?.stats >> (saveUmiBam(s) ? "${alignerDir(s)}/samtools_stats/" : null)
            s.transcriptome?.coord_sorted_samtools?.flagstat >> (saveUmiBam(s) ? "${alignerDir(s)}/samtools_stats/" : null)
            s.transcriptome?.coord_sorted_samtools?.idxstats >> (saveUmiBam(s) ? "${alignerDir(s)}/samtools_stats/" : null)
            s.transcriptome?.sorted_bam >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.transcriptome?.filtered_bam >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.prepare_for_rsem_log >> "${alignerDir(s)}/umitools/prepare_for_quantification_log/"
            s.transcriptomic_dedup_log >> "${alignerDir(s)}/${umiDedupToolDir(s)}/transcriptomic_dedup_log/"
            s.transcriptome?.dedup_bam >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.transcriptome?.sorted_bam_index >> (saveUmiBam(s) ? "${alignerDir(s)}/" : null)
            s.transcriptome?.stats >> "${alignerDir(s)}/samtools_stats/"
            s.transcriptome?.flagstat >> "${alignerDir(s)}/samtools_stats/"
            s.transcriptome?.idxstats >> "${alignerDir(s)}/samtools_stats/"
            s.transcriptome?.tsv?.edit_distance >> "${alignerDir(s)}/umitools/"
            s.transcriptome?.tsv?.per_umi >> "${alignerDir(s)}/umitools/"
            s.transcriptome?.tsv?.umi_per_position >> "${alignerDir(s)}/umitools/"
        }
    }

    markdup {   // MarkdupBam; anchor: metrics
        path { s ->
            s.metrics >> "${alignerDir(s)}/picard_metrics/"
            s.bam >> "${alignerDir(s)}/"
            s.cram >> "${alignerDir(s)}/"
            s.bai >> "${alignerDir(s)}/"
            s.samtools?.stats >> "${alignerDir(s)}/samtools_stats/"
            s.samtools?.flagstat >> "${alignerDir(s)}/samtools_stats/"
            s.samtools?.idxstats >> "${alignerDir(s)}/samtools_stats/"
        }
    }

    quant {   // RsemQuantSample | PseudoQuantSample (bam-salmon reuses the pseudo-alignment shape)
        path { s ->
            s.counts_gene >> "${alignerDir(s)}/"
            s.counts_transcript >> "${alignerDir(s)}/"
            s.stat >> "${alignerDir(s)}/"
            s.log >> "${alignerDir(s)}/log/"
            s.quant_dir >> "${alignerDir(s)}/"
        }
    }

    quant_merged {   // RsemQuantMerged | QuantMerged; never sample-prefixed except under --skip_quantification_merge
        path { r ->
            r.tpm_gene >> "${alignerDir(r)}/"
            r.counts_gene >> "${alignerDir(r)}/"
            r.lengths_gene >> "${alignerDir(r)}/"
            r.counts_gene_length_scaled >> (params.skip_quantification_merge ? null : "${alignerDir(r)}/")
            r.counts_gene_scaled >> "${alignerDir(r)}/"
            r.tpm_transcript >> "${alignerDir(r)}/"
            r.counts_transcript >> "${alignerDir(r)}/"
            r.lengths_transcript >> "${alignerDir(r)}/"
            r.tx2gene >> "${params.aligner}/"
            r.tx2gene_augmented >> "${alignerDir(r)}/"
            r.merged_gene_rds >> "${alignerDir(r)}/"
            r.merged_transcript_rds >> "${alignerDir(r)}/"
        }
    }

    quant_rsem_merge {   // record(id, rsem_merge: RsemMerge); rsem_merge is guaranteed non-null, ch_quant_rsem_merge filters out nulls before this target
        path { r ->
            r.rsem_merge.counts_gene >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge.tpm_gene >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge.counts_transcript >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge.tpm_transcript >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge.genes_long >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge.isoforms_long >> "${alignerDir(r)}/rsem_merge_counts/"
        }
    }

    quant_pseudo {   // PseudoQuantSample, pseudo-aligner
        path { s ->
            // s.log (kallisto only; always null for salmon) lives inside
            // quant_dir already and is not routed separately.
            s.quant_dir >> "${pseudoAlignerDir(s)}/"
        }
    }

    quant_merged_pseudo {   // QuantMerged, pseudo-aligner
        path { r ->
            r.tpm_gene >> "${pseudoAlignerDir(r)}/"
            r.counts_gene >> "${pseudoAlignerDir(r)}/"
            r.lengths_gene >> "${pseudoAlignerDir(r)}/"
            r.counts_gene_length_scaled >> (params.skip_quantification_merge ? null : "${pseudoAlignerDir(r)}/")
            r.counts_gene_scaled >> "${pseudoAlignerDir(r)}/"
            r.tpm_transcript >> "${pseudoAlignerDir(r)}/"
            r.counts_transcript >> "${pseudoAlignerDir(r)}/"
            r.lengths_transcript >> "${pseudoAlignerDir(r)}/"
            r.tx2gene >> "${params.pseudo_aligner}/"
            r.tx2gene_augmented >> "${pseudoAlignerDir(r)}/"
            r.merged_gene_rds >> "${pseudoAlignerDir(r)}/"
            r.merged_transcript_rds >> "${pseudoAlignerDir(r)}/"
        }
    }

    deseq2 {   // record(rdata, pca_vals, plots_pdf, sample_dists, size_factors, log); anchor: rdata; never sample-prefixed
        path { d ->
            d.rdata >> "${params.aligner}/deseq2_qc/"
            d.pca_vals >> "${params.aligner}/deseq2_qc/"
            d.plots_pdf >> "${params.aligner}/deseq2_qc/"
            d.sample_dists >> "${params.aligner}/deseq2_qc/"
            d.size_factors >> "${params.aligner}/deseq2_qc/"
            d.log >> "${params.aligner}/deseq2_qc/"
        }
    }

    deseq2_pseudo {   // same shape, pseudo-aligner; never sample-prefixed
        path { d ->
            d.rdata >> "${params.pseudo_aligner}/deseq2_qc/"
            d.pca_vals >> "${params.pseudo_aligner}/deseq2_qc/"
            d.plots_pdf >> "${params.pseudo_aligner}/deseq2_qc/"
            d.sample_dists >> "${params.pseudo_aligner}/deseq2_qc/"
            d.size_factors >> "${params.pseudo_aligner}/deseq2_qc/"
            d.log >> "${params.pseudo_aligner}/deseq2_qc/"
        }
    }

    bam_qc {   // BamQcRnaseq: preseq, featurecounts, biotype, qualimap, dupradar, rseqc
        enabled defineQcTools(params).size() > 0   // same check that decides whether any of these tools ran
        path { s ->
            s.preseq?.lc_extrap >> "${alignerDir(s)}/preseq/"
            s.preseq?.log >> "${alignerDir(s)}/preseq/log/"
            s.featurecounts?.counts >> "${alignerDir(s)}/featurecounts/"
            s.featurecounts?.summary >> "${alignerDir(s)}/featurecounts/"
            s.biotype?.tsv >> "${alignerDir(s)}/featurecounts/"
            s.biotype?.rrna >> "${alignerDir(s)}/featurecounts/"
            s.qualimap >> "${alignerDir(s)}/qualimap/"
            s.dupradar?.scatter2d >> "${alignerDir(s)}/dupradar/scatter_plot/"
            s.dupradar?.boxplot >> "${alignerDir(s)}/dupradar/box_plot/"
            s.dupradar?.hist >> "${alignerDir(s)}/dupradar/histogram/"
            s.dupradar?.dupmatrix >> "${alignerDir(s)}/dupradar/gene_data/"
            s.dupradar?.intercept_slope >> "${alignerDir(s)}/dupradar/intercepts_slope/"
            s.rseqc?.bamstat >> "${alignerDir(s)}/rseqc/bam_stat/"
            s.rseqc?.inferexperiment >> "${alignerDir(s)}/rseqc/infer_experiment/"
            s.rseqc?.junctionannotation?.pdf >> "${alignerDir(s)}/rseqc/junction_annotation/pdf/"
            s.rseqc?.junctionannotation?.events_pdf >> "${alignerDir(s)}/rseqc/junction_annotation/pdf/"
            s.rseqc?.junctionannotation?.bed >> "${alignerDir(s)}/rseqc/junction_annotation/bed/"
            s.rseqc?.junctionannotation?.interact_bed >> "${alignerDir(s)}/rseqc/junction_annotation/bed/"
            s.rseqc?.junctionannotation?.xls >> "${alignerDir(s)}/rseqc/junction_annotation/xls/"
            s.rseqc?.junctionannotation?.log >> "${alignerDir(s)}/rseqc/junction_annotation/log/"
            s.rseqc?.junctionannotation?.rscript >> "${alignerDir(s)}/rseqc/junction_annotation/rscript/"
            s.rseqc?.junctionsaturation?.pdf >> "${alignerDir(s)}/rseqc/junction_saturation/pdf/"
            s.rseqc?.junctionsaturation?.rscript >> "${alignerDir(s)}/rseqc/junction_saturation/rscript/"
            s.rseqc?.readdistribution >> "${alignerDir(s)}/rseqc/read_distribution/"
            s.rseqc?.readduplication?.pdf >> "${alignerDir(s)}/rseqc/read_duplication/pdf/"
            s.rseqc?.readduplication?.seq_xls >> "${alignerDir(s)}/rseqc/read_duplication/xls/"
            s.rseqc?.readduplication?.pos_xls >> "${alignerDir(s)}/rseqc/read_duplication/xls/"
            s.rseqc?.readduplication?.rscript >> "${alignerDir(s)}/rseqc/read_duplication/rscript/"
            s.rseqc?.innerdistance?.distance >> "${alignerDir(s)}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.freq >> "${alignerDir(s)}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.mean >> "${alignerDir(s)}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.pdf >> "${alignerDir(s)}/rseqc/inner_distance/pdf/"
            s.rseqc?.innerdistance?.rscript >> "${alignerDir(s)}/rseqc/inner_distance/rscript/"
            s.rseqc?.tin?.txt >> "${alignerDir(s)}/rseqc/tin/"
            s.rseqc?.tin?.xls >> "${alignerDir(s)}/rseqc/tin/"
        }
    }

    bam_qc_rustqc {   // record(id, file, target); --use_rustqc alternative, one entry per output file, target precomputed by rustqcTarget(); ch_bam_qc_rustqc_files filters out entries with a null target before this target
        path { s -> s.file >> s.target }
    }

    multiqc {   // MultiqcReport; anchor: report
        path { m ->
            m.report >> "${multiqcDir(m)}/"
            m.data >> "${multiqcDir(m)}/"
            m.plots >> "${multiqcDir(m)}/"
        }
    }

    pipeline_info {   // record(versions); anchor: versions
        path { p -> p.versions >> 'pipeline_info/' }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
