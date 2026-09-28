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
include { saveFile                   } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { outputEnabled              } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

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
        .flatMap { sample_id, r, _meta, runs, percent_mapped ->
            def transcriptome_bam = params.aligner == 'bowtie2_salmon' ? r.orig_bam : r.transcriptome_bam
            runs.collect { run ->
                record(
                    sample:            sample_id,
                    fastq_1:           run[0],
                    fastq_2:           run.size() > 1 ? run[1] : null,
                    strandedness:      r.meta.strandedness,
                    seq_platform:      r.meta.seq_platform ?: params.seq_platform,
                    seq_center:        r.meta.seq_center ?: params.seq_center,
                    genome_bam:        r.bam ? "${params.outdir}/${samplePrefix(r)}${params.aligner}/${r.bam.name}" : null,
                    percent_mapped:    percent_mapped,
                    transcriptome_bam: transcriptome_bam ? "${params.outdir}/${samplePrefix(r)}${params.aligner}/${transcriptome_bam.name}" : null
                )
            }
        }

    emit:
    trim_status         = RNASEQ.out.trim_status         // channel: [id, boolean]
    map_status          = RNASEQ.out.map_status          // channel: [id, boolean]
    strand_status       = RNASEQ.out.strand_status       // channel: [id, boolean]
    multiqc_report      = RNASEQ.out.multiqc_report      // channel: /path/to/multiqc_report.html
    genome              = ch_genome                      // channel: GenomeReferences fields + index: GenomeIndices + rrna_references
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
    genome           = NFCORE_RNASEQ.out.genome
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
            s.kraken2?.report                      >> "${samplePrefix(s)}${params.aligner}/contaminants/kraken2/kraken_reports/"
            s.kraken2?.classified_reads_fastq      >> (params.save_kraken_assignments ? "${samplePrefix(s)}${params.aligner}/contaminants/kraken2/kraken_reports/" : null)
            s.kraken2?.unclassified_reads_fastq    >> (params.save_kraken_assignments ? "${samplePrefix(s)}${params.aligner}/contaminants/kraken2/kraken_reports/" : null)
            s.kraken2?.classified_reads_assignment >> (params.save_kraken_unassigned ? "${samplePrefix(s)}${params.aligner}/contaminants/kraken2/kraken_reports/" : null)
            s.bracken?.abundance                   >> "${samplePrefix(s)}${params.aligner}/contaminants/bracken/"
            s.bracken?.report                      >> "${samplePrefix(s)}${params.aligner}/contaminants/bracken/"
            s.sylph?.profile                       >> "${samplePrefix(s)}${params.aligner}/contaminants/sylph/"
            s.sylphtax?.taxprof                    >> "${samplePrefix(s)}${params.aligner}/contaminants/sylph/"
        }
    }

    stringtie {   // record(id, meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo: StringtieAssembly?)
        enabled !params.skip_stringtie
        path { s ->
            s.transcript_gtf >> "${samplePrefix(s)}${params.aligner}/stringtie/"
            s.abundance >> "${samplePrefix(s)}${params.aligner}/stringtie/"
            s.coverage_gtf >> "${samplePrefix(s)}${params.aligner}/stringtie/"
            s.ballgown >> "${samplePrefix(s)}${params.aligner}/stringtie/"
            s.denovo?.transcript_gtf >> "${samplePrefix(s)}${params.aligner}/stringtie/"
        }
    }

    stringtie_merged {   // record(id, merged_gtf); cross-sample, only under --stringtie_ignore_gtf, never sample-prefixed
        enabled !params.skip_stringtie && params.stringtie_ignore_gtf
        path { s -> s.merged_gtf >> "${params.aligner}/stringtie/" }
    }

    bigwig {   // record(id, meta, combined, forward, reverse), each BigwigFiles { bigwig, bedgraph }; bedgraph stays unrouted
        enabled !params.skip_bigwig
        path { s ->
            s.combined?.bigwig >> "${samplePrefix(s)}${params.aligner}/bigwig/"
            s.forward?.bigwig  >> "${samplePrefix(s)}${params.aligner}/bigwig/"
            s.reverse?.bigwig  >> "${samplePrefix(s)}${params.aligner}/bigwig/"
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

    rrna_seqkit {   // record(bowtie2_index, seqkit_prefixed, seqkit_converted); ch_rrna_seqkit guarantees at least one FASTA list is non-empty
        path { r ->
            r.seqkit_prefixed >> 'seqkit/'
            r.seqkit_converted >> 'seqkit/'
        }
    }

    preprocessed {   // FastqQcTrimFilterSetstrandedness; no single anchor field survives every skip combination
        enabled outputEnabled('preprocessed')
        path { s ->
            s.fastqc?.raw_html >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.raw_zip >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.trim_html >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.trim_zip >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.filtered_html >> "${samplePrefix(s)}fastqc/filtered/"
            s.fastqc?.filtered_zip >> "${samplePrefix(s)}fastqc/filtered/"
            s.trim?.html >> "${samplePrefix(s)}${params.trimmer}/"
            s.trim?.log >> "${samplePrefix(s)}${params.trimmer == 'fastp' ? 'fastp/log/' : "${params.trimmer}/"}"
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
            // Only the single-end --un-gz FASTQs are published; paired-end reads rebuilt by SAMTOOLS_FASTQ_BOWTIE2 are not.
            s.reads >> (saveFile('non_ribo_reads', s) ? "${samplePrefix(s)}${params.ribo_removal_tool == 'bowtie2' ? 'bowtie2_rrna' : params.ribo_removal_tool}/" : null)
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
            s.samtools?.stats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.samtools?.flagstat >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.samtools?.idxstats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.bam >> (saveFile('align_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.bai >> (saveFile('align_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.orig_bam >> (params.save_align_intermeds ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome_bam >> (params.save_align_intermeds ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.unmapped >> (params.save_unaligned ? "${samplePrefix(s)}${params.aligner}/unmapped/" : null)
            s.star?.log_final >> "${samplePrefix(s)}${params.aligner}/log/"
            s.star?.log_out >> "${samplePrefix(s)}${params.aligner}/log/"
            s.star?.log_progress >> "${samplePrefix(s)}${params.aligner}/log/"
            s.star?.tab >> "${samplePrefix(s)}${params.aligner}/log/"
            s.hisat2?.summary >> "${samplePrefix(s)}${params.aligner}/log/"
            s.bowtie2?.log >> "${samplePrefix(s)}${params.aligner}/log/"
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
            s.genome?.stats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.genome?.flagstat >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.genome?.idxstats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.bam >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.bai >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.genomic_dedup_log >> "${samplePrefix(s)}${params.aligner}/${params.umi_dedup_tool}/genomic_dedup_log/"
            s.tsv?.edit_distance >> "${samplePrefix(s)}${params.aligner}/umitools/"
            s.tsv?.per_umi >> "${samplePrefix(s)}${params.aligner}/umitools/"
            s.tsv?.umi_per_position >> "${samplePrefix(s)}${params.aligner}/umitools/"
            s.transcriptome?.coord_sorted_bam >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome?.coord_sorted_bam_index >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome?.coord_sorted_samtools?.stats >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/samtools_stats/" : null)
            s.transcriptome?.coord_sorted_samtools?.flagstat >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/samtools_stats/" : null)
            s.transcriptome?.coord_sorted_samtools?.idxstats >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/samtools_stats/" : null)
            s.transcriptome?.sorted_bam >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome?.filtered_bam >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.prepare_for_rsem_log >> "${samplePrefix(s)}${params.aligner}/umitools/prepare_for_quantification_log/"
            s.transcriptomic_dedup_log >> "${samplePrefix(s)}${params.aligner}/${params.umi_dedup_tool}/transcriptomic_dedup_log/"
            s.transcriptome?.dedup_bam >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome?.sorted_bam_index >> (saveFile('umi_bam') ? "${samplePrefix(s)}${params.aligner}/" : null)
            s.transcriptome?.samtools?.stats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.transcriptome?.samtools?.flagstat >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.transcriptome?.samtools?.idxstats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.transcriptome?.tsv?.edit_distance >> "${samplePrefix(s)}${params.aligner}/umitools/"
            s.transcriptome?.tsv?.per_umi >> "${samplePrefix(s)}${params.aligner}/umitools/"
            s.transcriptome?.tsv?.umi_per_position >> "${samplePrefix(s)}${params.aligner}/umitools/"
        }
    }

    markdup {   // MarkdupBam; anchor: metrics
        path { s ->
            s.metrics >> "${samplePrefix(s)}${params.aligner}/picard_metrics/"
            s.bam >> "${samplePrefix(s)}${params.aligner}/"
            s.cram >> "${samplePrefix(s)}${params.aligner}/"
            s.bai >> "${samplePrefix(s)}${params.aligner}/"
            s.samtools?.stats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.samtools?.flagstat >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
            s.samtools?.idxstats >> "${samplePrefix(s)}${params.aligner}/samtools_stats/"
        }
    }

    quant {   // RsemQuantSample | PseudoQuantSample (bam-salmon reuses the pseudo-alignment shape)
        path { s ->
            s.counts_gene >> "${samplePrefix(s)}${params.aligner}/"
            s.counts_transcript >> "${samplePrefix(s)}${params.aligner}/"
            s.stat >> "${samplePrefix(s)}${params.aligner}/"
            s.log >> "${samplePrefix(s)}${params.aligner}/log/"
            s.quant_dir >> "${samplePrefix(s)}${params.aligner}/"
        }
    }

    quant_merged {   // QuantMerged; never sample-prefixed except under --skip_quantification_merge
        path { r ->
            r.tpm_gene >> "${samplePrefix(r)}${params.aligner}/"
            r.counts_gene >> "${samplePrefix(r)}${params.aligner}/"
            r.lengths_gene >> "${samplePrefix(r)}${params.aligner}/"
            r.counts_gene_length_scaled >> (params.skip_quantification_merge ? null : "${params.aligner}/")
            r.counts_gene_scaled >> "${samplePrefix(r)}${params.aligner}/"
            r.tpm_transcript >> "${samplePrefix(r)}${params.aligner}/"
            r.counts_transcript >> "${samplePrefix(r)}${params.aligner}/"
            r.lengths_transcript >> "${samplePrefix(r)}${params.aligner}/"
            r.tx2gene >> "${params.aligner}/"
            r.tx2gene_augmented >> "${samplePrefix(r)}${params.aligner}/"
            r.merged_gene_rds >> "${samplePrefix(r)}${params.aligner}/"
            r.merged_transcript_rds >> "${samplePrefix(r)}${params.aligner}/"
        }
    }

    // CUSTOM_RSEMMERGECOUNTS and tximport both write rsem.merged.* basenames, so they need separate targets (nextflow-io/nextflow#6617).
    quant_rsem_merge {   // record(id, rsem_merge: RsemMerge); empty channel unless --aligner star_rsem, rsem_merge is always non-null in every record it carries
        path { r ->
            r.rsem_merge.counts_gene >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
            r.rsem_merge.tpm_gene >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
            r.rsem_merge.counts_transcript >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
            r.rsem_merge.tpm_transcript >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
            r.rsem_merge.genes_long >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
            r.rsem_merge.isoforms_long >> "${samplePrefix(r)}${params.aligner}/rsem_merge_counts/"
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
            r.tpm_gene >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.counts_gene >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.lengths_gene >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.counts_gene_length_scaled >> (params.skip_quantification_merge ? null : "${params.pseudo_aligner}/")
            r.counts_gene_scaled >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.tpm_transcript >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.counts_transcript >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.lengths_transcript >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.tx2gene >> "${params.pseudo_aligner}/"
            r.tx2gene_augmented >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.merged_gene_rds >> "${samplePrefix(r)}${params.pseudo_aligner}/"
            r.merged_transcript_rds >> "${samplePrefix(r)}${params.pseudo_aligner}/"
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
            s.preseq?.lc_extrap >> "${samplePrefix(s)}${params.aligner}/preseq/"
            s.preseq?.log >> "${samplePrefix(s)}${params.aligner}/preseq/log/"
            s.featurecounts?.counts >> "${samplePrefix(s)}${params.aligner}/featurecounts/"
            s.featurecounts?.summary >> "${samplePrefix(s)}${params.aligner}/featurecounts/"
            s.biotype?.tsv >> "${samplePrefix(s)}${params.aligner}/featurecounts/"
            s.biotype?.rrna >> "${samplePrefix(s)}${params.aligner}/featurecounts/"
            s.qualimap >> "${samplePrefix(s)}${params.aligner}/qualimap/"
            s.dupradar?.scatter2d >> "${samplePrefix(s)}${params.aligner}/dupradar/scatter_plot/"
            s.dupradar?.boxplot >> "${samplePrefix(s)}${params.aligner}/dupradar/box_plot/"
            s.dupradar?.hist >> "${samplePrefix(s)}${params.aligner}/dupradar/histogram/"
            s.dupradar?.dupmatrix >> "${samplePrefix(s)}${params.aligner}/dupradar/gene_data/"
            s.dupradar?.intercept_slope >> "${samplePrefix(s)}${params.aligner}/dupradar/intercepts_slope/"
            s.rseqc?.bamstat >> "${samplePrefix(s)}${params.aligner}/rseqc/bam_stat/"
            s.rseqc?.inferexperiment >> "${samplePrefix(s)}${params.aligner}/rseqc/infer_experiment/"
            s.rseqc?.junctionannotation?.pdf >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/pdf/"
            s.rseqc?.junctionannotation?.events_pdf >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/pdf/"
            s.rseqc?.junctionannotation?.bed >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/bed/"
            s.rseqc?.junctionannotation?.interact_bed >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/bed/"
            s.rseqc?.junctionannotation?.xls >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/xls/"
            s.rseqc?.junctionannotation?.log >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/log/"
            s.rseqc?.junctionannotation?.rscript >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_annotation/rscript/"
            s.rseqc?.junctionsaturation?.pdf >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_saturation/pdf/"
            s.rseqc?.junctionsaturation?.rscript >> "${samplePrefix(s)}${params.aligner}/rseqc/junction_saturation/rscript/"
            s.rseqc?.readdistribution >> "${samplePrefix(s)}${params.aligner}/rseqc/read_distribution/"
            s.rseqc?.readduplication?.pdf >> "${samplePrefix(s)}${params.aligner}/rseqc/read_duplication/pdf/"
            s.rseqc?.readduplication?.seq_xls >> "${samplePrefix(s)}${params.aligner}/rseqc/read_duplication/xls/"
            s.rseqc?.readduplication?.pos_xls >> "${samplePrefix(s)}${params.aligner}/rseqc/read_duplication/xls/"
            s.rseqc?.readduplication?.rscript >> "${samplePrefix(s)}${params.aligner}/rseqc/read_duplication/rscript/"
            s.rseqc?.innerdistance?.distance >> "${samplePrefix(s)}${params.aligner}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.freq >> "${samplePrefix(s)}${params.aligner}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.mean >> "${samplePrefix(s)}${params.aligner}/rseqc/inner_distance/txt/"
            s.rseqc?.innerdistance?.pdf >> "${samplePrefix(s)}${params.aligner}/rseqc/inner_distance/pdf/"
            s.rseqc?.innerdistance?.rscript >> "${samplePrefix(s)}${params.aligner}/rseqc/inner_distance/rscript/"
            s.rseqc?.tin?.txt >> "${samplePrefix(s)}${params.aligner}/rseqc/tin/"
            s.rseqc?.tin?.xls >> "${samplePrefix(s)}${params.aligner}/rseqc/tin/"
        }
    }

    bam_qc_rustqc {   // record(id, file, target); --use_rustqc alternative, one entry per output file with a non-null target; the target is precomputed by rustqcTarget() in workflows/rnaseq/main.nf
        path { s -> s.file >> s.target }
    }

    multiqc {   // MultiqcReport; anchor: report
        path { m ->
            m.report >> "${samplePrefix(m)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"
            m.data >> "${samplePrefix(m)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"
            m.plots >> "${samplePrefix(m)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"
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
