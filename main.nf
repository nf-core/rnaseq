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

nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parameters that nextflow.config or the conf/ files read keep their defaults in
    nextflow.config, because the config is resolved before this block.
*/

params {
    // Input/output options
    input:                       Path? = null
    outdir:                      String? = null
    email:                       String? = null

    // Reference genome options
    fasta:                       String? = getGenomeAttribute('fasta')
    gtf:                         String? = getGenomeAttribute('gtf')
    gff:                         String? = getGenomeAttribute('gff')
    gene_bed:                    String? = getGenomeAttribute('bed12')
    transcript_fasta:            String? = getGenomeAttribute('transcript_fasta')
    additional_fasta:            String? = getGenomeAttribute('additional_fasta')
    splicesites:                 String? = null
    star_index:                  String? = getGenomeAttribute('star')
    hisat2_index:                String? = getGenomeAttribute('hisat2')
    rsem_index:                  String? = getGenomeAttribute('rsem')
    salmon_index:                String? = getGenomeAttribute('salmon')
    kallisto_index:              String? = getGenomeAttribute('kallisto')
    bowtie2_index:               String? = getGenomeAttribute('bowtie2')
    hisat2_build_memory:         String = '200.GB'
    gencode:                     Boolean = false
    prokaryotic:                 Boolean = false
    gffread_transcript_fasta:    Boolean = false
    gtf_extra_attributes:        String = 'gene_name'
    gtf_group_features:          String = 'gene_id'
    featurecounts_group_type:    String = 'gene_biotype'

    // Read trimming options
    trimmer:                     String = 'trimgalore'
    min_trimmed_reads:           Integer = 10000

    // Read filtering options
    bbsplit_fasta_list:          String? = null
    bbsplit_index:               String? = getGenomeAttribute('bbsplit')
    sortmerna_index:             String? = getGenomeAttribute('sortmerna')
    remove_ribo_rna:             Boolean = false
    ribo_removal_tool:           String = 'sortmerna'
    bowtie2_rrna_index:          String? = null
    ribo_database_manifest:      String = "${projectDir}/assets/rrna-db-defaults.txt"

    // UMI options
    with_umi:                    Boolean = false
    umi_dedup_tool:              String = 'umitools'
    umi_discard_read:            Integer = 0
    umitools_dedup_stats:        Boolean = false
    umitools_dedup_primary_only: Boolean = false

    // Alignment options
    aligner:                     String = 'star_salmon'
    use_sentieon_star:           Boolean = false
    use_parabricks_star:         Boolean = false
    pseudo_aligner:              String? = null
    bam_csi_index:               Boolean = false
    star_ignore_sjdbgtf:         Boolean = false
    min_mapped_reads:            Float = 5.0
    seq_center:                  String? = null
    seq_platform:                String? = null
    stringtie_ignore_gtf:        Boolean = false
    kallisto_quant_fraglen:      Integer = 200
    kallisto_quant_fraglen_sd:   Integer = 200
    stranded_threshold:          Float = 0.8
    unstranded_threshold:        Float = 0.1

    // Optional outputs
    save_merged_fastq:           Boolean = false
    save_umi_intermeds:          Boolean = false
    save_non_ribo_reads:         Boolean = false
    save_bbsplit_reads:          Boolean = false
    save_reference:              Boolean = false
    save_trimmed:                Boolean = false
    save_align_intermeds:        Boolean = false
    save_unaligned:              Boolean = false
    save_kraken_assignments:     Boolean = false
    save_kraken_unassigned:      Boolean = false

    // Quality Control
    rseqc_modules:               String = 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication'
    contaminant_screening:       String? = null
    contaminant_screening_input: String = 'unmapped'
    kraken_db:                   String? = null
    sylph_db:                    Path? = null
    sylph_taxonomy:              Path? = null

    // Process skipping options
    skip_gtf_filter:             Boolean = false
    skip_bbsplit:                Boolean = true
    skip_umi_extract:            Boolean = false
    skip_linting:                Boolean = false
    skip_trimming:               Boolean = false
    skip_alignment:              Boolean = false
    skip_pseudo_alignment:       Boolean = false
    skip_quantification_merge:   Boolean = false
    skip_markduplicates:         Boolean = false
    skip_bigwig:                 Boolean = false
    skip_stringtie:              Boolean = false
    skip_fastqc:                 Boolean = false
    use_rustqc:                  Boolean = false
    skip_preseq:                 Boolean = true
    skip_dupradar:               Boolean = false
    skip_qualimap:               Boolean = false
    skip_rseqc:                  Boolean = false
    skip_biotype_qc:             Boolean = false
    skip_deseq2_qc:              Boolean = false
    skip_multiqc:                Boolean = false
    skip_qc:                     Boolean = false

    // Generic options
    version:                     Boolean = false
    email_on_fail:               String? = null
    plaintext_email:             Boolean = false
    monochrome_logs:             Boolean = false
    multiqc_config:              Path? = null
    multiqc_logo:                Path? = null
    multiqc_methods_description: Path? = null
    validate_params:             Boolean = true
    help:                        String? = null
    help_full:                   Boolean = false
    show_hidden:                 Boolean = false
}

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
include { getGenomeAttribute         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy          } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { anySampleAutoStrandedness  } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

include { AlignedSample; Contaminants; StringtieSample; BigwigSample; QuantSample; RsemMergeSample; Deseq2Qc; RustqcResult; LintFile; PipelineInfo; SamplesheetRow } from './subworkflows/local/types'
include { GenomeArtifact                                   } from './subworkflows/local/utils_nfcore_rnaseq_pipeline/types'
include { FastqQcTrimFilterSetstrandedness; RrnaReferences } from './subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/types'
include { UmiDedupBam                                      } from './subworkflows/nf-core/bam_dedup_umi/types'
include { MarkdupBam                                       } from './subworkflows/nf-core/bam_markduplicates_picard/types'
include { BamQcRnaseq                                      } from './subworkflows/nf-core/bam_qc_rnaseq/types'
include { QuantMerged                                      } from './subworkflows/nf-core/quant_tximport_summarizedexperiment/types'
include { StringtieMerged                                  } from './subworkflows/nf-core/bam_stringtie_merge/types'
include { MultiqcReport                                    } from './subworkflows/local/multiqc_rnaseq/types'

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
    def references_ribo_tool: String? = params.remove_ribo_rna && !(params.ribo_removal_tool == "bowtie2" && params.bowtie2_rrna_index) ? params.ribo_removal_tool as String : null
    def indices_ribo_tool: String? = params.remove_ribo_rna ? params.ribo_removal_tool as String : null

    references = PREPARE_GENOME_REFERENCES (
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
        references_ribo_tool,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.contaminant_screening,
        params.prokaryotic
    )

    //
    // SUBWORKFLOW: Build or load aligner / pseudo-aligner / filtering indices
    //
    indices = PREPARE_GENOME_INDICES (
        references.fasta_fai,
        references.gtf,
        references.transcript_fasta,
        references.rrna_fastas,
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
        indices_ribo_tool,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.use_parabricks_star,
        isStarIndexLegacy() ? true : false,
        params.hisat2_build_memory,
        anySampleAutoStrandedness()
    )

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        references.fasta_fai.map { _meta, _fasta, fai -> fai != null ? checkMaxContigSize(fai) : null }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    qc_tools = defineQcTools(params)

    results = RNASEQ (
        ch_samplesheet,
        references.fasta_fai,
        references.gtf,
        references.chrom_sizes,
        references.gene_bed,
        references.transcript_fasta,
        indices.star_index,
        indices.rsem_index,
        indices.hisat2_index,
        indices.bowtie2_index,
        indices.salmon_index,
        indices.kallisto_index,
        indices.bbsplit_index,
        references.rrna_fastas,
        indices.sortmerna_index,
        indices.bowtie2_rrna_index,
        indices.splicesites,
        references.kraken_db,
        qc_tools
    )

    // Matches current behavior: no task-workdir guard, so a user-supplied
    // --bowtie2_rrna_index is republished here exactly as it is today.
    ch_rrna_bowtie2_index = results.rrna_references
        .map { r -> r.bowtie2_index }
        .flatMap { index -> index == null ? [] : [index] }

    // Same-basename fields split out per stage to avoid a >> rename-key collision (nextflow-io/nextflow#6617).
    ch_lint_raw     = results.preprocessed.map { r -> record(id: r.id, file: r.lint?.raw) }.filter { s -> s.file != null }
    ch_lint_trimmed = results.preprocessed.map { r -> record(id: r.id, file: r.lint?.trimmed) }.filter { s -> s.file != null }
    ch_lint_bbsplit = results.preprocessed.map { r -> record(id: r.id, file: r.lint?.bbsplit) }.filter { s -> s.file != null }
    ch_lint_ribo    = results.preprocessed.map { r -> record(id: r.id, file: r.lint?.ribo) }.filter { s -> s.file != null }

    // The prepared rRNA FASTAs publish whether or not --save_reference is set, so they cannot ride the genome target.
    ch_rrna_seqkit = results.rrna_references.flatMap { r -> (r.seqkit_prefixed ?: []).isEmpty() && (r.seqkit_converted ?: []).isEmpty() ? [] : [r] }

    // samplesheet_with_bams.csv rows: one per sequencing run of a sample, with the aligned record's
    // meta so an inferred strandedness replaces 'auto'. genome_bam is the
    // coordinate-sorted BAM; bowtie2_salmon aligns to the transcriptome, so its unsorted bowtie2 BAM is the transcriptome_bam.
    // The BAMs are published by the `aligned` output, so the row carries their published location
    // as a string; routing the Paths through `samplesheet` too would copy each file twice.
    ch_runs           = results.reads.map { meta, runs -> record(id: meta.id, runs: runs) }
    ch_percent_mapped = results.percent_mapped.map { id, percent_mapped -> record(id: id, percent_mapped: percent_mapped) }

    ch_samplesheet_rows = results.aligned
        .join(ch_runs, by: 'id')
        .join(ch_percent_mapped, by: 'id')
        .flatMap { r ->
            r.runs.collect { run ->
                def transcriptome_bam = params.aligner == 'bowtie2_salmon' ? r.orig_bam[0] : r.transcriptome_bam
                record(
                    sample:            r.id,
                    fastq_1:           run[0],
                    fastq_2:           run.size() > 1 ? run[1] : null,
                    strandedness:      r.meta.strandedness,
                    seq_platform:      r.meta.seq_platform ?: params.seq_platform,
                    seq_center:        r.meta.seq_center ?: params.seq_center,
                    genome_bam:        r.bam != null ? "${params.outdir}/${alignedDir(r.id)}${r.bam.name}" : null,
                    percent_mapped:    r.percent_mapped,
                    transcriptome_bam: transcriptome_bam != null ? "${params.outdir}/${alignedDir(r.id)}${transcriptome_bam.name}" : null
                )
            }
        }

    emit:
    trim_status:          Channel<Tuple<String, Boolean>>           = results.trim_status
    map_status:           Channel<Tuple<String, Boolean>>           = results.map_status
    strand_status:        Channel<Tuple<String, Boolean>>           = results.strand_status
    multiqc_report:       Channel<Path>                             = results.multiqc_report
    genome_references:    Channel<GenomeArtifact>                   = references.references
    genome_intermediates: Channel<GenomeArtifact>                   = references.intermediates
    genome_indices:       Channel<GenomeArtifact>                   = indices.indices
    rrna_bowtie2_index:   Channel<Path>                             = ch_rrna_bowtie2_index
    rrna_seqkit:          Channel<RrnaReferences>                   = ch_rrna_seqkit

    // Stage result records, keyed on id
    preprocessed:         Channel<FastqQcTrimFilterSetstrandedness> = results.preprocessed
    lint_raw:             Channel<LintFile>                         = ch_lint_raw
    lint_trimmed:         Channel<LintFile>                         = ch_lint_trimmed
    lint_bbsplit:         Channel<LintFile>                         = ch_lint_bbsplit
    lint_ribo:            Channel<LintFile>                         = ch_lint_ribo
    aligned:              Channel<AlignedSample>                    = results.aligned
    umi_dedup:            Channel<UmiDedupBam>                      = results.umi_dedup
    markdup:              Channel<MarkdupBam>                       = results.markdup
    bam_qc:               Channel<BamQcRnaseq>                      = results.bam_qc
    bam_qc_rustqc:        Channel<RustqcResult>                     = results.bam_qc_rustqc
    samplesheet:          Channel<SamplesheetRow>                   = ch_samplesheet_rows
    quant:                Channel<QuantSample>                      = results.quant
    quant_merged:         Channel<QuantMerged>                      = results.quant_merged
    quant_rsem_merge:     Channel<RsemMergeSample>                  = results.quant_rsem_merge
    quant_pseudo:         Channel<QuantSample>                      = results.quant_pseudo
    quant_merged_pseudo:  Channel<QuantMerged>                      = results.quant_merged_pseudo
    contaminants:         Channel<Contaminants>               = results.contaminants
    stringtie:            Channel<StringtieSample>                  = results.stringtie
    bigwig:               Channel<BigwigSample>                     = results.bigwig

    // Run-level result records
    stringtie_merged:     Channel<StringtieMerged>                  = results.stringtie_merged
    deseq2:               Channel<Deseq2Qc>                    = results.deseq2
    deseq2_pseudo:        Channel<Deseq2Qc>                    = results.deseq2_pseudo
    multiqc:              Channel<MultiqcReport>                    = results.multiqc
    pipeline_info:        Channel<PipelineInfo>                     = results.pipeline_info
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
    results = NFCORE_RNASEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        results.multiqc_report,
        results.trim_status,
        results.map_status,
        results.strand_status
    )

    publish:
    contaminants         = results.contaminants
    stringtie            = results.stringtie
    stringtie_merged     = results.stringtie_merged
    bigwig               = results.bigwig
    genome_references    = results.genome_references
    genome_intermediates = results.genome_intermediates
    genome_indices       = results.genome_indices
    rrna_bowtie2_index   = results.rrna_bowtie2_index
    rrna_seqkit          = results.rrna_seqkit
    preprocessed         = results.preprocessed
    lint_raw             = results.lint_raw
    lint_trimmed         = results.lint_trimmed
    lint_bbsplit         = results.lint_bbsplit
    lint_ribo            = results.lint_ribo
    aligned              = results.aligned
    umi_dedup            = results.umi_dedup
    markdup              = results.markdup
    samplesheet          = results.samplesheet
    quant                = results.quant
    quant_merged         = results.quant_merged
    quant_rsem_merge     = results.quant_rsem_merge
    quant_pseudo         = results.quant_pseudo
    quant_merged_pseudo  = results.quant_merged_pseudo
    deseq2               = results.deseq2
    deseq2_pseudo        = results.deseq2_pseudo
    bam_qc               = results.bam_qc
    bam_qc_rustqc        = results.bam_qc_rustqc
    multiqc              = results.multiqc
    pipeline_info        = results.pipeline_info
}

// Per-sample directory prefix for per-record outputs under --skip_quantification_merge.
// Run-level records (e.g. a cross-sample merged file) are never sample-prefixed.
def samplePrefix(id: String) -> String { params.skip_quantification_merge ? "${id}/" : '' }

// Directory the aligned BAM/BAI/transcriptome BAM publish to. Also used to rebuild the
// samplesheet's genome_bam/transcriptome_bam paths in the entry workflow, so the two can't drift.
def alignedDir(id: String) -> String { "${samplePrefix(id)}${params.aligner}/" }

output {
    contaminants: Channel<Contaminants> {   // record(id, meta, kraken2, bracken, sylph, sylphtax); exactly one tool branch is populated per run
        enabled !params.skip_qc && params.contaminant_screening
        path { s ->
            def dir     = "${alignedDir(s.id)}contaminants/"
            def kraken  = "${dir}kraken2/kraken_reports/"
            def bracken = "${dir}bracken/"
            def sylph   = "${dir}sylph/"
            def saveAssignments = params.save_kraken_assignments ? kraken : null
            [
                (s.kraken2?.report):                         kraken,
                (s.kraken2?.classified_reads_fastq):         saveAssignments,
                (s.kraken2?.unclassified_reads_fastq):       saveAssignments,
                (s.kraken2?.classified_reads_assignment):    params.save_kraken_unassigned ? kraken : null,
                ([s.bracken?.abundance, s.bracken?.report]): bracken,
                ([s.sylph?.profile, s.sylphtax?.taxprof]):   sylph,
            ]
        }
    }

    stringtie: Channel<StringtieSample> {   // record(id, meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo: StringtieAssembly?)
        enabled !params.skip_stringtie
        path { s ->
            def dir = "${alignedDir(s.id)}stringtie/"
            [
                ([s.transcript_gtf, s.abundance, s.coverage_gtf, s.denovo?.transcript_gtf]): dir,
                (s.ballgown):                                                                dir,
            ]
        }
    }

    stringtie_merged: Channel<StringtieMerged> {   // record(id, merged_gtf); cross-sample, only under --stringtie_ignore_gtf, never sample-prefixed
        enabled !params.skip_stringtie && params.stringtie_ignore_gtf
        path { s -> s.merged_gtf >> "${params.aligner}/stringtie/" }
    }

    bigwig: Channel<BigwigSample> {   // record(id, meta, combined, forward, reverse), each BigwigFiles { bigwig, bedgraph }; bedgraph stays unrouted
        enabled !params.skip_bigwig
        path { s ->
            [([s.combined?.bigwig, s.forward?.bigwig, s.reverse?.bigwig]): "${alignedDir(s.id)}bigwig/"]
        }
    }

    genome_references: Channel<GenomeArtifact> {   // GenomeArtifact stream, one record per top-level reference file; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> (r.kind == 'kraken_db' ? 'genome/index/' : 'genome/') }
    }

    genome_intermediates: Channel<GenomeArtifact> {   // GenomeArtifact stream, one record per superseded/incidental reference file; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> 'genome/' }
    }

    genome_indices: Channel<GenomeArtifact> {   // GenomeArtifact stream, one record per index/log; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> (r.kind == 'sortmerna' ? 'genome/sortmerna/' : 'genome/index/') }
    }

    rrna_bowtie2_index: Channel<Path> {   // path(bowtie2_rrna/index/); never sample-prefixed
        enabled params.save_reference
        path { p -> p >> 'bowtie2_rrna/index/' }
    }

    rrna_seqkit: Channel<RrnaReferences> {   // record(bowtie2_index, seqkit_prefixed, seqkit_converted); ch_rrna_seqkit guarantees at least one FASTA list is non-empty
        path { r ->
            r.seqkit_prefixed >> 'seqkit/'
            r.seqkit_converted >> 'seqkit/'
        }
    }

    preprocessed: Channel<FastqQcTrimFilterSetstrandedness> {   // FastqQcTrimFilterSetstrandedness; no single anchor field survives every skip combination
        path { s ->
            def sp       = samplePrefix(s.id)
            def fastqc   = "${sp}fastqc/"
            def trimDir  = "${sp}${params.trimmer}/"
            def umitools = "${sp}umitools/"
            def bbsplit  = "${sp}bbsplit/"
            def riboDir  = "${sp}${params.ribo_removal_tool == 'bowtie2' ? 'bowtie2_rrna' : params.ribo_removal_tool}/"
            def isFastp  = params.trimmer == 'fastp'
            def saveTrimmed = params.save_trimmed ? trimDir : null
            def saveBbsplit = params.save_bbsplit_reads ? bbsplit : null
            // Only the single-end --un-gz FASTQs are published; paired-end reads rebuilt by SAMTOOLS_FASTQ_BOWTIE2 are not.
            def saveNonRibo = params.remove_ribo_rna && params.save_non_ribo_reads && (params.ribo_removal_tool != 'bowtie2' || s.meta.single_end)
            [
                // `reads` can be the same files as `reads_cat` or `reads_trimmed`. A later entry for the same files
                // replaces an earlier one even when its target is null, so `reads` must stay first.
                (s.reads):                       saveNonRibo ? riboDir : null,
                (s.reads_cat):                   params.save_merged_fastq ? "${sp}fastq/" : null,
                (s.reads_trimmed):               params.skip_trimming ? null : saveTrimmed,
                (s.fastqc?.raw_html):            "${fastqc}raw/",
                (s.fastqc?.raw_zip):             "${fastqc}raw/",
                (s.fastqc?.trim_html):           "${fastqc}trim/",
                (s.fastqc?.trim_zip):            "${fastqc}trim/",
                (s.fastqc?.filtered_html):       "${fastqc}filtered/",
                (s.fastqc?.filtered_zip):        "${fastqc}filtered/",
                (s.trim?.html):                  trimDir,
                (s.trim?.log):                   isFastp ? "${trimDir}log/" : trimDir,
                (s.trim?.json):                  isFastp ? trimDir : null,
                (s.trim?.unpaired):              saveTrimmed,
                (s.trim?.reads_fail):            saveTrimmed,
                (s.trim?.reads_merged):          saveTrimmed,
                (s.umi?.log):                    umitools,
                (s.umi?.reads):                  params.save_umi_intermeds ? umitools : null,
                (s.bbsplit?.stats):              bbsplit,
                (s.bbsplit?.primary_reads):      saveBbsplit,
                (s.bbsplit?.other_genome_reads): saveBbsplit,
                (s.rrna?.sortmerna_log):         "${sp}sortmerna/",
                ([s.rrna?.ribodetector_log, s.rrna?.seqkit_stats]): "${sp}ribodetector/",
                (s.rrna?.bowtie2_log):           "${sp}bowtie2_rrna/",
            ]
        }
    }

    // One target per FQ_LINT stage: same basename, avoids the >> collision (nextflow-io/nextflow#6617).
    lint_raw: Channel<LintFile> {   // record(id, file); file is guaranteed non-null, ch_lint_raw filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s.id)}fq_lint/raw/" }
    }

    lint_trimmed: Channel<LintFile> {   // record(id, file); file is guaranteed non-null, ch_lint_trimmed filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s.id)}fq_lint/trimmed/" }
    }

    lint_bbsplit: Channel<LintFile> {   // record(id, file); file is guaranteed non-null, ch_lint_bbsplit filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s.id)}fq_lint/bbsplit/" }
    }

    lint_ribo: Channel<LintFile> {   // record(id, file); file is guaranteed non-null, ch_lint_ribo filters out nulls before this target
        path { s -> s.file >> "${samplePrefix(s.id)}fq_lint/${params.ribo_removal_tool ?: 'sortmerna'}/" }
    }

    aligned: Channel<AlignedSample> {   // STAR, Bowtie2 or HISAT2 result; anchor: samtools.stats
        path { s ->
            def dir         = alignedDir(s.id)
            def logDir      = "${dir}log/"
            def bamDir      = (params.save_align_intermeds || params.skip_markduplicates) ? dir : null
            def intermedDir = params.save_align_intermeds ? dir : null
            [
                ([s.samtools?.stats, s.samtools?.flagstat, s.samtools?.idxstats]): "${dir}samtools_stats/",
                ([s.bam, s.bai]):                                                   bamDir,
                (s.orig_bam):                                                       intermedDir,
                (s.transcriptome_bam):                                              intermedDir,
                (s.unmapped):                                                       params.save_unaligned ? "${dir}unmapped/" : null,
                ([s.star?.log_final, s.star?.log_out, s.star?.log_progress, s.hisat2?.summary, s.bowtie2?.log]): logDir,
                (s.star?.tab):                                                      logDir,
            ]
        }
    }

    samplesheet: Channel<SamplesheetRow> {   // record(sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam); field order is the CSV column order
        enabled params.save_align_intermeds && !params.skip_alignment
        index {
            path 'samplesheets/samplesheet_with_bams.csv'
            header true
        }
    }

    umi_dedup: Channel<UmiDedupBam> {   // UmiDedupBam; anchor: genome.stats
        path { s ->
            def dir      = alignedDir(s.id)
            def stats    = "${dir}samtools_stats/"
            def umitools = "${dir}umitools/"
            def tool     = "${dir}${params.umi_dedup_tool}/"
            def saveBam  = params.save_align_intermeds || params.save_umi_intermeds
            def t        = s.transcriptome
            [
                ([s.genome?.stats, s.genome?.flagstat, s.genome?.idxstats,
                  t?.samtools?.stats, t?.samtools?.flagstat, t?.samtools?.idxstats]): stats,
                ([s.bam, s.bai, t?.coord_sorted_bam, t?.coord_sorted_bam_index,
                  t?.sorted_bam, t?.sorted_bam_index, t?.filtered_bam, t?.dedup_bam]): saveBam ? dir : null,
                ([t?.coord_sorted_samtools?.stats, t?.coord_sorted_samtools?.flagstat,
                  t?.coord_sorted_samtools?.idxstats]): saveBam ? stats : null,
                ([s.tsv?.edit_distance, s.tsv?.per_umi, s.tsv?.umi_per_position,
                  t?.tsv?.edit_distance, t?.tsv?.per_umi, t?.tsv?.umi_per_position]): umitools,
                (s.prepare_for_rsem_log):        "${umitools}prepare_for_quantification_log/",
                (s.genomic_dedup_log):           "${tool}genomic_dedup_log/",
                (s.transcriptomic_dedup_log):    "${tool}transcriptomic_dedup_log/",
            ]
        }
    }

    markdup: Channel<MarkdupBam> {   // MarkdupBam; anchor: metrics
        path { s ->
            def dir = alignedDir(s.id)
            [
                ([s.bam, s.cram, s.bai]):                                           dir,
                (s.metrics):                                                        "${dir}picard_metrics/",
                ([s.samtools?.stats, s.samtools?.flagstat, s.samtools?.idxstats]): "${dir}samtools_stats/",
            ]
        }
    }

    quant: Channel<QuantSample> {   // RSEM or Salmon-on-BAM result (bam-salmon reuses the pseudo-alignment shape)
        path { s ->
            def dir = alignedDir(s.id)
            [
                ([s.counts_gene, s.counts_transcript, s.stat, s.quant_dir]): dir,
                (s.log):                                                     "${dir}log/",
            ]
        }
    }

    quant_merged: Channel<QuantMerged> {   // QuantMerged; never sample-prefixed except under --skip_quantification_merge
        path { r ->
            def top = "${params.aligner}/"
            def dir = "${samplePrefix(r.id)}${top}"
            [
                ([r.tpm_gene, r.counts_gene, r.lengths_gene, r.counts_gene_scaled,
                  r.tpm_transcript, r.counts_transcript, r.lengths_transcript,
                  r.tx2gene_augmented, r.merged_gene_rds, r.merged_transcript_rds]): dir,
                (r.counts_gene_length_scaled): params.skip_quantification_merge ? null : top,
                (r.tx2gene):                   top,
            ]
        }
    }

    // CUSTOM_RSEMMERGECOUNTS and tximport both write rsem.merged.* basenames, so they need separate targets (nextflow-io/nextflow#6617).
    quant_rsem_merge: Channel<RsemMergeSample> {   // record(id, rsem_merge: RsemMerge); empty channel unless --aligner star_rsem, rsem_merge is always non-null in every record it carries
        path { r ->
            def m = r.rsem_merge
            [
                ([m.counts_gene, m.tpm_gene, m.counts_transcript, m.tpm_transcript, m.genes_long, m.isoforms_long]): "${alignedDir(r.id)}rsem_merge_counts/",
            ]
        }
    }

    quant_pseudo: Channel<QuantSample> {   // pseudo-aligner result
        path { s ->
            // s.log (kallisto only; always null for salmon) lives inside
            // quant_dir already and is not routed separately.
            s.quant_dir >> "${samplePrefix(s.id)}${params.pseudo_aligner}/"
        }
    }

    quant_merged_pseudo: Channel<QuantMerged> {   // QuantMerged, pseudo-aligner
        path { r ->
            def top = "${params.pseudo_aligner}/"
            def dir = "${samplePrefix(r.id)}${top}"
            [
                ([r.tpm_gene, r.counts_gene, r.lengths_gene, r.counts_gene_scaled,
                  r.tpm_transcript, r.counts_transcript, r.lengths_transcript,
                  r.tx2gene_augmented, r.merged_gene_rds, r.merged_transcript_rds]): dir,
                (r.counts_gene_length_scaled): params.skip_quantification_merge ? null : top,
                (r.tx2gene):                   top,
            ]
        }
    }

    deseq2: Channel<Deseq2Qc> {   // record(rdata, pca_vals, plots_pdf, sample_dists, size_factors, log); anchor: rdata; never sample-prefixed
        path { d ->
            [([d.rdata, d.pca_vals, d.plots_pdf, d.sample_dists, d.size_factors, d.log]): "${params.aligner}/deseq2_qc/"]
        }
    }

    deseq2_pseudo: Channel<Deseq2Qc> {   // same shape, pseudo-aligner; never sample-prefixed
        path { d ->
            [([d.rdata, d.pca_vals, d.plots_pdf, d.sample_dists, d.size_factors, d.log]): "${params.pseudo_aligner}/deseq2_qc/"]
        }
    }

    bam_qc: Channel<BamQcRnaseq> {   // BamQcRnaseq: preseq, featurecounts, biotype, qualimap, dupradar, rseqc
        enabled defineQcTools(params).size() > 0   // same check that decides whether any of these tools ran
        path { s ->
            def dir      = alignedDir(s.id)
            def dupradar = "${dir}dupradar/"
            def rseqc    = "${dir}rseqc/"
            def annDir   = "${rseqc}junction_annotation/"
            def satDir   = "${rseqc}junction_saturation/"
            def rdupDir  = "${rseqc}read_duplication/"
            def innerDir = "${rseqc}inner_distance/"
            def rq       = s.rseqc
            [
                (s.preseq?.lc_extrap):          "${dir}preseq/",
                (s.preseq?.log):                "${dir}preseq/log/",
                ([s.featurecounts?.counts, s.featurecounts?.summary, s.biotype?.tsv, s.biotype?.rrna]): "${dir}featurecounts/",
                (s.qualimap):                   "${dir}qualimap/",
                (s.dupradar?.scatter2d):        "${dupradar}scatter_plot/",
                (s.dupradar?.boxplot):          "${dupradar}box_plot/",
                (s.dupradar?.hist):             "${dupradar}histogram/",
                (s.dupradar?.dupmatrix):        "${dupradar}gene_data/",
                (s.dupradar?.intercept_slope):  "${dupradar}intercepts_slope/",
                (rq?.bamstat):                  "${rseqc}bam_stat/",
                (rq?.inferexperiment):          "${rseqc}infer_experiment/",
                ([rq?.junctionannotation?.pdf, rq?.junctionannotation?.events_pdf]):  "${annDir}pdf/",
                ([rq?.junctionannotation?.bed, rq?.junctionannotation?.interact_bed]): "${annDir}bed/",
                (rq?.junctionannotation?.xls):     "${annDir}xls/",
                (rq?.junctionannotation?.log):     "${annDir}log/",
                (rq?.junctionannotation?.rscript): "${annDir}rscript/",
                (rq?.junctionsaturation?.pdf):     "${satDir}pdf/",
                (rq?.junctionsaturation?.rscript): "${satDir}rscript/",
                (rq?.readdistribution):            "${rseqc}read_distribution/",
                (rq?.readduplication?.pdf):        "${rdupDir}pdf/",
                ([rq?.readduplication?.seq_xls, rq?.readduplication?.pos_xls]): "${rdupDir}xls/",
                (rq?.readduplication?.rscript):    "${rdupDir}rscript/",
                ([rq?.innerdistance?.distance, rq?.innerdistance?.freq, rq?.innerdistance?.mean]): "${innerDir}txt/",
                (rq?.innerdistance?.pdf):          "${innerDir}pdf/",
                (rq?.innerdistance?.rscript):      "${innerDir}rscript/",
                ([rq?.tin?.txt, rq?.tin?.xls]):    "${rseqc}tin/",
            ]
        }
    }

    bam_qc_rustqc: Channel<RustqcResult> {   // RUSTQC record (meta, samtools, preseq, dupradar, featurecounts, biotype, rseqc, qualimap); --use_rustqc alternative to bam_qc, never sample-prefixed
        path { s ->
            def dir      = "${params.aligner}/rustqc/"
            def dupradar = "${dir}dupradar/"
            def fcDir    = "${dir}featurecounts/"
            def rseqc    = "${dir}rseqc/"
            def annDir   = "${rseqc}junction_annotation/"
            def satDir   = "${rseqc}junction_saturation/"
            def rdupDir  = "${rseqc}read_duplication/"
            def innerDir = "${rseqc}inner_distance/"
            def rq       = s.rseqc
            [
                ([s.samtools?.stats, s.samtools?.flagstat, s.samtools?.idxstats]): "${dir}samtools_stats/",
                (s.preseq?.lc_extrap):          "${dir}preseq/",
                (s.dupradar?.scatter2d):        "${dupradar}scatter_plot/",
                (s.dupradar?.boxplot):          "${dupradar}box_plot/",
                (s.dupradar?.hist):             "${dupradar}histogram/",
                (s.dupradar?.dupmatrix):        "${dupradar}gene_data/",
                (s.dupradar?.intercept_slope):  "${dupradar}intercepts_slope/",
                (s.dupradar?.multiqc):          dupradar,
                ([s.featurecounts?.counts, s.biotype?.tsv, s.biotype?.mqc, s.biotype?.rrna]): fcDir,
                // Published under the name featureCounts gives its non-biotype summary
                (s.featurecounts?.summary):     "${fcDir}${s.featurecounts?.summary?.name?.replace('.biotype.tsv.summary', '.tsv.summary')}",
                (s.qualimap):                   "${dir}qualimap/",
                (rq?.bamstat):                  "${rseqc}bam_stat/",
                (rq?.inferexperiment):          "${rseqc}infer_experiment/",
                (rq?.readdistribution):         "${rseqc}read_distribution/",
                ([rq?.tin?.txt, rq?.tin?.xls]): "${rseqc}tin/",
                ([rq?.junctionannotation?.bed, rq?.junctionannotation?.interact_bed]): "${annDir}bed/",
                (rq?.junctionannotation?.xls):     "${annDir}xls/",
                (rq?.junctionannotation?.log):     "${annDir}log/",
                (rq?.junctionannotation?.plot):    "${annDir}plot/",
                (rq?.junctionannotation?.rscript): "${annDir}rscript/",
                (rq?.junctionsaturation?.summary): "${satDir}txt/",
                (rq?.junctionsaturation?.plot):    "${satDir}plot/",
                (rq?.junctionsaturation?.rscript): "${satDir}rscript/",
                ([rq?.readduplication?.seq_xls, rq?.readduplication?.pos_xls]): "${rdupDir}xls/",
                (rq?.readduplication?.plot):       "${rdupDir}plot/",
                (rq?.readduplication?.rscript):    "${rdupDir}rscript/",
                ([rq?.innerdistance?.distance, rq?.innerdistance?.freq, rq?.innerdistance?.mean, rq?.innerdistance?.summary]): "${innerDir}txt/",
                (rq?.innerdistance?.plot):         "${innerDir}plot/",
                (rq?.innerdistance?.rscript):      "${innerDir}rscript/",
            ]
        }
    }

    multiqc: Channel<MultiqcReport> {   // MultiqcReport; anchor: report
        path { m ->
            [([m.report, m.data, m.plots]): "${samplePrefix(m.id)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"]
        }
    }

    pipeline_info: Channel<PipelineInfo> {   // record(versions); anchor: versions
        path { p -> p.versions >> 'pipeline_info/' }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
