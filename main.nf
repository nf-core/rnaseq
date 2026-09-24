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

    emit:
    trim_status         = RNASEQ.out.trim_status         // channel: [id, boolean]
    map_status          = RNASEQ.out.map_status          // channel: [id, boolean]
    strand_status       = RNASEQ.out.strand_status       // channel: [id, boolean]
    multiqc_report      = RNASEQ.out.multiqc_report      // channel: /path/to/multiqc_report.html
    genome              = ch_genome                      // channel: GenomeReferences fields + index: GenomeIndices + rrna_references

    // Stage result records, keyed on id
    preprocessed        = RNASEQ.out.preprocessed        // channel: FastqQcTrimFilterSetstrandedness
    aligned             = RNASEQ.out.aligned             // channel: StarAligned | Bowtie2Aligned | Hisat2Aligned
    umi_dedup           = RNASEQ.out.umi_dedup           // channel: UmiDedupBam
    markdup             = RNASEQ.out.markdup             // channel: MarkdupBam
    bam_qc              = RNASEQ.out.bam_qc              // channel: BamQcRnaseq
    bam_qc_rustqc       = RNASEQ.out.bam_qc_rustqc       // channel: record(id, meta, samtools, dupradar, featurecounts, preseq, rseqc, qualimap)
    quant               = RNASEQ.out.quant               // channel: RsemQuantSample | PseudoQuantSample, alignment-based quantifier
    quant_merged        = RNASEQ.out.quant_merged        // channel: RsemQuantMerged | QuantMerged, alignment-based quantifier
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
    // Routed one area at a time; unlisted channels still use conf/modules/*.config publishDir.
    contaminants     = NFCORE_RNASEQ.out.contaminants
    stringtie        = NFCORE_RNASEQ.out.stringtie
    stringtie_merged = NFCORE_RNASEQ.out.stringtie_merged
    bigwig           = NFCORE_RNASEQ.out.bigwig
    genome           = NFCORE_RNASEQ.out.genome
    preprocessed     = NFCORE_RNASEQ.out.preprocessed
    aligned          = NFCORE_RNASEQ.out.aligned
    umi_dedup        = NFCORE_RNASEQ.out.umi_dedup
    markdup          = NFCORE_RNASEQ.out.markdup
    quant               = NFCORE_RNASEQ.out.quant
    quant_merged        = NFCORE_RNASEQ.out.quant_merged
    quant_pseudo        = NFCORE_RNASEQ.out.quant_pseudo
    quant_merged_pseudo = NFCORE_RNASEQ.out.quant_merged_pseudo
    deseq2              = NFCORE_RNASEQ.out.deseq2
    deseq2_pseudo       = NFCORE_RNASEQ.out.deseq2_pseudo
    bam_qc              = NFCORE_RNASEQ.out.bam_qc
    multiqc             = NFCORE_RNASEQ.out.multiqc
    pipeline_info       = NFCORE_RNASEQ.out.pipeline_info
}

// Run-level records (e.g. a cross-sample merged file) are never sample-prefixed.
def samplePrefix(r) { params.skip_quantification_merge ? "${r.id}/" : '' }
def alignerDir(r)    { "${samplePrefix(r)}${params.aligner}" }
def trimLogDir(s)    { params.trimmer == 'fastp' ? "${samplePrefix(s)}${params.trimmer}/log" : "${samplePrefix(s)}${params.trimmer}" }
def saveAlignBam(s)  { params.save_align_intermeds || params.skip_markduplicates }
def saveUmiBam(s)    { params.save_align_intermeds || params.save_umi_intermeds }
def umiDedupToolDir(s) { params.umi_dedup_tool == 'umicollapse' ? 'umicollapse' : 'umitools' }
def pseudoAlignerDir(r) { "${samplePrefix(r)}${params.pseudo_aligner}" }
def multiqcDir(m) {
    def suffix = params.skip_alignment ? '' : "/${params.aligner}"
    m.id == 'multiqc_report' ? "multiqc${suffix}" : "${m.id}/multiqc${suffix}"
}

// The dir the surviving reads would have published to under the mechanism
// that last touched them (rRNA removal, then BBSplit), or null if neither
// ran - trimming alone never publishes preprocessed.reads on its own.
def readsLastStageDir(s) {
    if (s.rrna?.sortmerna_log != null)    { return params.save_non_ribo_reads ? "${samplePrefix(s)}sortmerna" : null }
    if (s.rrna?.ribodetector_log != null) { return params.save_non_ribo_reads ? "${samplePrefix(s)}ribodetector" : null }
    if (s.rrna?.bowtie2_log != null)      { return params.save_non_ribo_reads ? "${samplePrefix(s)}bowtie2_rrna" : null }
    if (s.bbsplit != null)                { return params.save_bbsplit_reads ? "${samplePrefix(s)}bbsplit" : null }
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
            g.index?.sortmerna >> 'genome/sortmerna/'
            g.index?.bowtie2_rrna >> 'genome/index/'
            g.rrna_references?.bowtie2_index >> 'bowtie2_rrna/index/'
        }
    }

    preprocessed {   // FastqQcTrimFilterSetstrandedness; anchor: lint.raw (unless --skip_linting) or fastqc.raw_zip
        path { s ->
            s.fastqc?.raw_html >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.raw_zip >> "${samplePrefix(s)}fastqc/raw/"
            s.fastqc?.trim_html >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.trim_zip >> "${samplePrefix(s)}fastqc/trim/"
            s.fastqc?.filtered_html >> "${samplePrefix(s)}fastqc/filtered/"
            s.fastqc?.filtered_zip >> "${samplePrefix(s)}fastqc/filtered/"
            s.trim?.html >> "${samplePrefix(s)}${params.trimmer}/"
            s.trim?.log >> "${trimLogDir(s)}/"
            s.trim?.json >> "${samplePrefix(s)}${params.trimmer}/"
            s.trim?.unpaired >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.trim?.reads_fail >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.trim?.reads_merged >> (params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.reads_trimmed >> (!params.skip_trimming && params.save_trimmed ? "${samplePrefix(s)}${params.trimmer}/" : null)
            s.umi?.log >> "${samplePrefix(s)}umitools/"
            s.umi?.reads >> (params.save_umi_intermeds ? "${samplePrefix(s)}umitools/" : null)
            s.bbsplit?.stats >> "${samplePrefix(s)}bbsplit/"
            s.bbsplit?.other_genome_reads >> (params.save_bbsplit_reads ? "${samplePrefix(s)}bbsplit/" : null)
            s.lint?.raw >> "${samplePrefix(s)}fq_lint/raw/"
            s.lint?.trimmed >> "${samplePrefix(s)}fq_lint/trimmed/"
            s.lint?.bbsplit >> "${samplePrefix(s)}fq_lint/bbsplit/"
            s.lint?.ribo >> "${samplePrefix(s)}fq_lint/${params.ribo_removal_tool ?: 'sortmerna'}/"
            s.rrna?.sortmerna_log >> "${samplePrefix(s)}sortmerna/"
            s.rrna?.ribodetector_log >> "${samplePrefix(s)}ribodetector/"
            s.rrna?.seqkit_stats >> "${samplePrefix(s)}ribodetector/"
            s.rrna?.bowtie2_log >> "${samplePrefix(s)}bowtie2_rrna/"
            s.reads_cat >> (params.save_merged_fastq ? "${samplePrefix(s)}fastq/" : null)
            s.reads >> readsLastStageDir(s)
        }
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
            r.tx2gene >> "${alignerDir(r)}/"
            r.tx2gene_augmented >> "${alignerDir(r)}/"
            r.merged_gene_rds >> "${alignerDir(r)}/"
            r.merged_transcript_rds >> "${alignerDir(r)}/"
            r.rsem_merge?.counts_gene >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge?.tpm_gene >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge?.counts_transcript >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge?.tpm_transcript >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge?.genes_long >> "${alignerDir(r)}/rsem_merge_counts/"
            r.rsem_merge?.isoforms_long >> "${alignerDir(r)}/rsem_merge_counts/"
        }
    }

    quant_pseudo {   // PseudoQuantSample, pseudo-aligner
        path { s ->
            s.quant_dir >> "${pseudoAlignerDir(s)}/"
            s.log >> "${pseudoAlignerDir(s)}/"
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
            r.tx2gene >> "${pseudoAlignerDir(r)}/"
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

    // bam_qc_rustqc (the --use_rustqc alternative) is not routed here yet: RustQC's
    // own saveAs closure renames some files (e.g. the biotype summary), which >>
    // cannot do - it only chooses a directory, never a filename. Needs a decision
    // (module-level rename vs. accepting a path/name change) before it's routed.
    bam_qc {   // BamQcRnaseq: preseq, featurecounts, biotype, qualimap, dupradar, rseqc
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

    multiqc {   // MultiqcReport; anchor: report
        path { m ->
            m.report >> "${multiqcDir(m)}/"
            m.data >> "${multiqcDir(m)}/"
            m.plots >> "${multiqcDir(m)}/"
        }
    }

    pipeline_info {   // record(versions); anchor: versions. Needs VM verification:
                       // collectFile() without storeDir must produce a task-output
                       // path for >> to route (nextflow-io/nextflow#7667).
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
