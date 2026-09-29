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

include { samplesheetToList          } from 'plugin/nf-schema'
include { RNASEQ                     } from './workflows/rnaseq'
include { PREPARE_GENOME_REFERENCES  } from './subworkflows/local/prepare_genome_references'
include { PREPARE_GENOME_INDICES     } from './subworkflows/local/prepare_genome_indices'
include { PIPELINE_INITIALISATION    } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION        } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { defineQcTools              } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy          } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { samplePrefix               } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { alignedDir                 } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

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
        params.hisat2_build_memory,
        anySampleAutoStrandedness()
    )

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

    // Matches current behavior: no task-workdir guard, so a user-supplied
    // --bowtie2_rrna_index is republished here exactly as it is today.
    ch_rrna_bowtie2_index = RNASEQ.out.rrna_references.map { r -> r.bowtie2_index }

    // Same-basename fields split out per stage to avoid a >> rename-key collision (nextflow-io/nextflow#6617).
    ch_lint_raw     = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.raw) }.filter { s -> s.file != null }
    ch_lint_trimmed = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.trimmed) }.filter { s -> s.file != null }
    ch_lint_bbsplit = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.bbsplit) }.filter { s -> s.file != null }
    ch_lint_ribo    = RNASEQ.out.preprocessed.map { r -> record(id: r.id, file: r.lint?.ribo) }.filter { s -> s.file != null }

    // The prepared rRNA FASTAs publish whether or not --save_reference is set, so they cannot ride the genome target.
    ch_rrna_seqkit = RNASEQ.out.rrna_references.filter { r -> r.seqkit_prefixed || r.seqkit_converted }

    // samplesheet_with_bams.csv rows: one per sequencing run of a sample, with the aligned record's
    // meta so an inferred strandedness replaces 'auto'. genome_bam is the
    // coordinate-sorted BAM; bowtie2_salmon aligns to the transcriptome, so its unsorted bowtie2 BAM is the transcriptome_bam.
    // The BAMs are published by the `aligned` output, so the row carries their published location
    // as a string; routing the Paths through `samplesheet` too would copy each file twice.
    ch_samplesheet_rows = RNASEQ.out.aligned
        .map { r -> [r.id, r] }
        .join(RNASEQ.out.reads.map { meta, runs -> [meta.id, meta, runs] })
        .join(RNASEQ.out.percent_mapped)
        .transpose(by: 3) // one row per sequencing run: replicates sample_id/r/_meta/percent_mapped across each entry in runs
        .map { sample_id, r, _meta, run, percent_mapped ->
            def transcriptome_bam = params.aligner == 'bowtie2_salmon' ? r.orig_bam : r.transcriptome_bam
            record(
                sample:            sample_id,
                fastq_1:           run[0],
                fastq_2:           run.size() > 1 ? run[1] : null,
                strandedness:      r.meta.strandedness,
                seq_platform:      r.meta.seq_platform ?: params.seq_platform,
                seq_center:        r.meta.seq_center ?: params.seq_center,
                genome_bam:        r.bam ? "${params.outdir}/${alignedDir(r)}${r.bam.name}" : null,
                percent_mapped:    percent_mapped,
                transcriptome_bam: transcriptome_bam ? "${params.outdir}/${alignedDir(r)}${transcriptome_bam.name}" : null
            )
        }

    emit:
    trim_status         = RNASEQ.out.trim_status         // channel: [id, boolean]
    map_status          = RNASEQ.out.map_status          // channel: [id, boolean]
    strand_status       = RNASEQ.out.strand_status       // channel: [id, boolean]
    multiqc_report      = RNASEQ.out.multiqc_report      // channel: /path/to/multiqc_report.html
    genome_references   = PREPARE_GENOME_REFERENCES.out.references     // channel: GenomeArtifact, one record per top-level reference file actually built or supplied
    genome_intermediates = PREPARE_GENOME_REFERENCES.out.intermediates // channel: GenomeArtifact, one record per superseded/incidental reference file
    genome_indices      = PREPARE_GENOME_INDICES.out.indices           // channel: GenomeArtifact, one record per index/log actually built or supplied
    rrna_bowtie2_index  = ch_rrna_bowtie2_index                        // channel: path(bowtie2_rrna/index/), only when the bowtie2 rRNA index is built
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
    bam_qc_rustqc       = RNASEQ.out.bam_qc_rustqc       // channel: record(id, file, target), one entry per RustQC output file
    samplesheet         = ch_samplesheet_rows            // channel: record(sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam), one entry per sequencing run
    quant               = RNASEQ.out.quant               // channel: RsemQuantSample | PseudoQuantSample, alignment-based quantifier
    quant_merged        = RNASEQ.out.quant_merged        // channel: QuantMerged, alignment-based quantifier
    quant_rsem_merge    = RNASEQ.out.quant_rsem_merge    // channel: record(id, rsem_merge: RsemMerge), CUSTOM_RSEMMERGECOUNTS outputs; empty unless --aligner star_rsem
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
    genome_references    = NFCORE_RNASEQ.out.genome_references
    genome_intermediates = NFCORE_RNASEQ.out.genome_intermediates
    genome_indices       = NFCORE_RNASEQ.out.genome_indices
    rrna_bowtie2_index   = NFCORE_RNASEQ.out.rrna_bowtie2_index
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

output {
    contaminants {   // record(id, meta, kraken2, bracken, sylph, sylphtax); exactly one tool branch is populated per run
        enabled !params.skip_qc && params.contaminant_screening
        path { s ->
            def dir     = "${alignedDir(s)}contaminants/"
            def kraken  = "${dir}kraken2/kraken_reports/"
            def bracken = "${dir}bracken/"
            def sylph   = "${dir}sylph/"
            [
                (s.kraken2?.report):                         kraken,
                (s.kraken2?.classified_reads_fastq):         params.save_kraken_assignments ? kraken : null,
                (s.kraken2?.unclassified_reads_fastq):       params.save_kraken_assignments ? kraken : null,
                (s.kraken2?.classified_reads_assignment):    params.save_kraken_unassigned ? kraken : null,
                ([s.bracken?.abundance, s.bracken?.report]): bracken,
                ([s.sylph?.profile, s.sylphtax?.taxprof]):   sylph,
            ]
        }
    }

    stringtie {   // record(id, meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo: StringtieAssembly?)
        enabled !params.skip_stringtie
        path { s ->
            def dir = "${alignedDir(s)}stringtie/"
            [
                ([s.transcript_gtf, s.abundance, s.coverage_gtf, s.denovo?.transcript_gtf]): dir,
                (s.ballgown):                                                                dir,
            ]
        }
    }

    stringtie_merged {   // record(id, merged_gtf); cross-sample, only under --stringtie_ignore_gtf, never sample-prefixed
        enabled !params.skip_stringtie && params.stringtie_ignore_gtf
        path { s -> s.merged_gtf >> "${params.aligner}/stringtie/" }
    }

    bigwig {   // record(id, meta, combined, forward, reverse), each BigwigFiles { bigwig, bedgraph }; bedgraph stays unrouted
        enabled !params.skip_bigwig
        path { s ->
            [([s.combined?.bigwig, s.forward?.bigwig, s.reverse?.bigwig]): "${alignedDir(s)}bigwig/"]
        }
    }

    genome_references {   // GenomeArtifact stream, one record per top-level reference file; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> (r.kind == 'kraken_db' ? 'genome/index/' : 'genome/') }
    }

    genome_intermediates {   // GenomeArtifact stream, one record per superseded/incidental reference file; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> 'genome/' }
    }

    genome_indices {   // GenomeArtifact stream, one record per index/log; never sample-prefixed
        enabled params.save_reference
        path { r -> r.file >> (r.kind == 'sortmerna' ? 'genome/sortmerna/' : 'genome/index/') }
    }

    rrna_bowtie2_index {   // path(bowtie2_rrna/index/); never sample-prefixed
        enabled params.save_reference
        path { p -> p >> 'bowtie2_rrna/index/' }
    }

    rrna_seqkit {   // record(bowtie2_index, seqkit_prefixed, seqkit_converted); ch_rrna_seqkit guarantees at least one FASTA list is non-empty
        path { r ->
            r.seqkit_prefixed >> 'seqkit/'
            r.seqkit_converted >> 'seqkit/'
        }
    }

    preprocessed {   // FastqQcTrimFilterSetstrandedness; no single anchor field survives every skip combination
        path { s ->
            def sp      = samplePrefix(s)
            def fastqc  = "${sp}fastqc/"
            def trimDir = "${sp}${params.trimmer}/"
            def riboDir = "${sp}${params.ribo_removal_tool == 'bowtie2' ? 'bowtie2_rrna' : params.ribo_removal_tool}/"
            // Only the single-end --un-gz FASTQs are published; paired-end reads rebuilt by SAMTOOLS_FASTQ_BOWTIE2 are not.
            def saveNonRibo = params.remove_ribo_rna && params.save_non_ribo_reads && (params.ribo_removal_tool != 'bowtie2' || s.meta.single_end)
            [
                // `reads` can be the same files as `reads_cat` or `reads_trimmed`. A later entry for the same files
                // replaces an earlier one even when its target is null, so `reads` must stay first.
                (s.reads):                       saveNonRibo ? riboDir : null,
                (s.reads_cat):                   params.save_merged_fastq ? "${sp}fastq/" : null,
                (s.reads_trimmed):               !params.skip_trimming && params.save_trimmed ? trimDir : null,
                (s.fastqc?.raw_html):            "${fastqc}raw/",
                (s.fastqc?.raw_zip):             "${fastqc}raw/",
                (s.fastqc?.trim_html):           "${fastqc}trim/",
                (s.fastqc?.trim_zip):            "${fastqc}trim/",
                (s.fastqc?.filtered_html):       "${fastqc}filtered/",
                (s.fastqc?.filtered_zip):        "${fastqc}filtered/",
                (s.trim?.html):                  trimDir,
                (s.trim?.log):                   params.trimmer == 'fastp' ? "${trimDir}log/" : trimDir,
                (s.trim?.json):                  params.trimmer == 'fastp' ? trimDir : null,
                (s.trim?.unpaired):              params.save_trimmed ? trimDir : null,
                (s.trim?.reads_fail):            params.save_trimmed ? trimDir : null,
                (s.trim?.reads_merged):          params.save_trimmed ? trimDir : null,
                (s.umi?.log):                    "${sp}umitools/",
                (s.umi?.reads):                  params.save_umi_intermeds ? "${sp}umitools/" : null,
                (s.bbsplit?.stats):              "${sp}bbsplit/",
                (s.bbsplit?.primary_reads):      params.save_bbsplit_reads ? "${sp}bbsplit/" : null,
                (s.bbsplit?.other_genome_reads): params.save_bbsplit_reads ? "${sp}bbsplit/" : null,
                (s.rrna?.sortmerna_log):         "${sp}sortmerna/",
                ([s.rrna?.ribodetector_log, s.rrna?.seqkit_stats]): "${sp}ribodetector/",
                (s.rrna?.bowtie2_log):           "${sp}bowtie2_rrna/",
            ]
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
            def dir     = alignedDir(s)
            def saveBam = params.save_align_intermeds || params.skip_markduplicates
            [
                ([s.samtools?.stats, s.samtools?.flagstat, s.samtools?.idxstats]): "${dir}samtools_stats/",
                ([s.bam, s.bai]):                                                   saveBam ? dir : null,
                (s.orig_bam):                                                       params.save_align_intermeds ? dir : null,
                (s.transcriptome_bam):                                              params.save_align_intermeds ? dir : null,
                (s.unmapped):                                                       params.save_unaligned ? "${dir}unmapped/" : null,
                ([s.star?.log_final, s.star?.log_out, s.star?.log_progress, s.hisat2?.summary, s.bowtie2?.log]): "${dir}log/",
                (s.star?.tab):                                                      "${dir}log/",
            ]
        }
    }

    samplesheet {   // record(sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam); field order is the CSV column order
        enabled params.save_align_intermeds && !params.skip_alignment
        index {
            path 'samplesheets/samplesheet_with_bams.csv'
            header true
        }
    }

    umi_dedup {   // UmiDedupBam; anchor: genome.stats
        path { s ->
            def dir      = alignedDir(s)
            def stats    = "${dir}samtools_stats/"
            def umitools = "${dir}umitools/"
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
                (s.genomic_dedup_log):           "${dir}${params.umi_dedup_tool}/genomic_dedup_log/",
                (s.transcriptomic_dedup_log):    "${dir}${params.umi_dedup_tool}/transcriptomic_dedup_log/",
            ]
        }
    }

    markdup {   // MarkdupBam; anchor: metrics
        path { s ->
            def dir = alignedDir(s)
            [
                ([s.bam, s.cram, s.bai]):                                           dir,
                (s.metrics):                                                        "${dir}picard_metrics/",
                ([s.samtools?.stats, s.samtools?.flagstat, s.samtools?.idxstats]): "${dir}samtools_stats/",
            ]
        }
    }

    quant {   // RsemQuantSample | PseudoQuantSample (bam-salmon reuses the pseudo-alignment shape)
        path { s ->
            def dir = alignedDir(s)
            [
                ([s.counts_gene, s.counts_transcript, s.stat, s.quant_dir]): dir,
                (s.log):                                                     "${dir}log/",
            ]
        }
    }

    quant_merged {   // QuantMerged; never sample-prefixed except under --skip_quantification_merge
        path { r ->
            def top = "${params.aligner}/"
            def dir = "${samplePrefix(r)}${top}"
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
    quant_rsem_merge {   // record(id, rsem_merge: RsemMerge); empty channel unless --aligner star_rsem, rsem_merge is always non-null in every record it carries
        path { r ->
            def m = r.rsem_merge
            [
                ([m.counts_gene, m.tpm_gene, m.counts_transcript, m.tpm_transcript, m.genes_long, m.isoforms_long]): "${alignedDir(r)}rsem_merge_counts/",
            ]
        }
    }

    quant_pseudo {   // PseudoQuantSample, pseudo-aligner
        path { s ->
            // s.log (kallisto only; always null for salmon) lives inside
            // quant_dir already and is not routed separately.
            s.quant_dir >> "${samplePrefix(s)}${params.pseudo_aligner}/"
        }
    }

    quant_merged_pseudo {   // QuantMerged, pseudo-aligner
        path { r ->
            def top = "${params.pseudo_aligner}/"
            def dir = "${samplePrefix(r)}${top}"
            [
                ([r.tpm_gene, r.counts_gene, r.lengths_gene, r.counts_gene_scaled,
                  r.tpm_transcript, r.counts_transcript, r.lengths_transcript,
                  r.tx2gene_augmented, r.merged_gene_rds, r.merged_transcript_rds]): dir,
                (r.counts_gene_length_scaled): params.skip_quantification_merge ? null : top,
                (r.tx2gene):                   top,
            ]
        }
    }

    deseq2 {   // record(rdata, pca_vals, plots_pdf, sample_dists, size_factors, log); anchor: rdata; never sample-prefixed
        path { d ->
            [([d.rdata, d.pca_vals, d.plots_pdf, d.sample_dists, d.size_factors, d.log]): "${params.aligner}/deseq2_qc/"]
        }
    }

    deseq2_pseudo {   // same shape, pseudo-aligner; never sample-prefixed
        path { d ->
            [([d.rdata, d.pca_vals, d.plots_pdf, d.sample_dists, d.size_factors, d.log]): "${params.pseudo_aligner}/deseq2_qc/"]
        }
    }

    bam_qc {   // BamQcRnaseq: preseq, featurecounts, biotype, qualimap, dupradar, rseqc
        enabled defineQcTools(params).size() > 0   // same check that decides whether any of these tools ran
        path { s ->
            def dir   = alignedDir(s)
            def dup   = "${dir}dupradar/"
            def rseqc = "${dir}rseqc/"
            def rq    = s.rseqc
            [
                (s.preseq?.lc_extrap):          "${dir}preseq/",
                (s.preseq?.log):                "${dir}preseq/log/",
                ([s.featurecounts?.counts, s.featurecounts?.summary, s.biotype?.tsv, s.biotype?.rrna]): "${dir}featurecounts/",
                (s.qualimap):                   "${dir}qualimap/",
                (s.dupradar?.scatter2d):        "${dup}scatter_plot/",
                (s.dupradar?.boxplot):          "${dup}box_plot/",
                (s.dupradar?.hist):             "${dup}histogram/",
                (s.dupradar?.dupmatrix):        "${dup}gene_data/",
                (s.dupradar?.intercept_slope):  "${dup}intercepts_slope/",
                (rq?.bamstat):                  "${rseqc}bam_stat/",
                (rq?.inferexperiment):          "${rseqc}infer_experiment/",
                ([rq?.junctionannotation?.pdf, rq?.junctionannotation?.events_pdf]):  "${rseqc}junction_annotation/pdf/",
                ([rq?.junctionannotation?.bed, rq?.junctionannotation?.interact_bed]): "${rseqc}junction_annotation/bed/",
                (rq?.junctionannotation?.xls):     "${rseqc}junction_annotation/xls/",
                (rq?.junctionannotation?.log):     "${rseqc}junction_annotation/log/",
                (rq?.junctionannotation?.rscript): "${rseqc}junction_annotation/rscript/",
                (rq?.junctionsaturation?.pdf):     "${rseqc}junction_saturation/pdf/",
                (rq?.junctionsaturation?.rscript): "${rseqc}junction_saturation/rscript/",
                (rq?.readdistribution):            "${rseqc}read_distribution/",
                (rq?.readduplication?.pdf):        "${rseqc}read_duplication/pdf/",
                ([rq?.readduplication?.seq_xls, rq?.readduplication?.pos_xls]): "${rseqc}read_duplication/xls/",
                (rq?.readduplication?.rscript):    "${rseqc}read_duplication/rscript/",
                ([rq?.innerdistance?.distance, rq?.innerdistance?.freq, rq?.innerdistance?.mean]): "${rseqc}inner_distance/txt/",
                (rq?.innerdistance?.pdf):          "${rseqc}inner_distance/pdf/",
                (rq?.innerdistance?.rscript):      "${rseqc}inner_distance/rscript/",
                ([rq?.tin?.txt, rq?.tin?.xls]):    "${rseqc}tin/",
            ]
        }
    }

    bam_qc_rustqc {   // record(id, file, target); --use_rustqc alternative, one entry per output file with a non-null target; the target is precomputed by rustqcTarget() in workflows/rnaseq/main.nf
        path { s -> s.file >> s.target }
    }

    multiqc {   // MultiqcReport; anchor: report
        path { m ->
            [([m.report, m.data, m.plots]): "${samplePrefix(m)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"]
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

//
// Check whether any sample declares strandedness 'auto'
//

def anySampleAutoStrandedness() {
    samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
        .any { meta, _fastq_1, _fastq_2, _genome_bam, _transcriptome_bam -> meta.strandedness == 'auto' }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
