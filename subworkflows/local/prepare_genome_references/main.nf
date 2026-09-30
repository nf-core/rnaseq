nextflow.enable.types = true

//
// Uncompress and prepare reference genome files (FASTA / GTF / BED / transcript FASTA / chrom.sizes / rRNA / Kraken DB)
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_RRNA_FASTAS      } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_KRAKEN_DB          } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { GFFREAD as GFFREAD_TRANSCRIPTS    } from '../../../modules/nf-core/gffread'
include { GFFREAD as GFFREAD_GENE_BED       } from '../../../modules/nf-core/gffread'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'
include { SENTIEON_RSEMPREPAREREFERENCE as SENTIEON_MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/sentieon/rsempreparereference'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { EAUTILS_GTF2BED                      } from '../../../modules/nf-core/ea-utils/gtf2bed'
include { CUSTOM_GTFFILTER                     } from '../../../modules/nf-core/custom/gtffilter'

include { taskOutputOrNull                     } from '../utils_nfcore_rnaseq_pipeline'
include { GenomeArtifact                       } from '../utils_nfcore_rnaseq_pipeline/types'

workflow PREPARE_GENOME_REFERENCES {

    take:
    fasta: String?                    // file: /path/to/genome.fasta (optional!)
    gtf: String?                      // file: /path/to/genome.gtf
    gff: String?                      // file: /path/to/genome.gff
    additional_fasta: String?         // file: /path/to/additional.fasta
    transcript_fasta: String?         // file: /path/to/transcript.fasta
    gene_bed: String?                 // file: /path/to/gene.bed
    sortmerna_fasta_list: String?     // file: /path/to/sortmerna_fasta_list.txt
    kraken_db: String?                // path: /path/to/kraken2/db/ or .tar.gz archive
    gencode: Boolean                  // whether the genome is from GENCODE
    gffread_transcript_fasta: Boolean // use gffread instead of RSEM for transcript FASTA extraction
    featurecounts_group_type: String  // The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner: String?                  // Specifies the alignment algorithm to use
    pseudo_aligner: String?           // Specifies the pseudo aligner to use
    skip_gtf_filter: Boolean          // Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    ribo_removal_tool: String?        // Tool for rRNA removal - 'sortmerna', 'ribodetector', or 'bowtie2' (null if skip)
    skip_alignment: Boolean           // Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment: Boolean    // Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star: Boolean        // whether to use sentieon STAR version
    contaminant_screening: String?    // contaminant screening tool ('kraken2', 'kraken2_bracken', 'sylph', or null)
    prokaryotic: Boolean              // whether the genome is prokaryotic (CDS-only annotation - use gffread --bed for gene BED since ea-utils/gtf2bed only handles exon features)

    main:
    // Absent artifacts are Values holding null so every stream keeps one type; steps that
    // dereference an artifact are guarded by the plain booleans below.
    ch_no_path = channel.value(null as Path)
    def has_gtf = (gtf || gff) ? true : false
    def fasta_provided = (fasta ? true : false)

    //---------------------------
    // 1) Uncompress GTF or GFF -> GTF
    //---------------------------
    ch_gff_uncompressed = ch_no_path
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF(record(id: 'gtf', meta: [:], archive: file(gtf, checkIfExists: true))).map { r -> r.gunzip }
        } else {
            ch_gtf = channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff_gunzipped    = GUNZIP_GFF(record(id: 'gff', meta: [:], archive: file(gff, checkIfExists: true)))
            ch_gff_uncompressed = ch_gff_gunzipped.map { r -> r.gunzip }
            ch_gff              = ch_gff_gunzipped.map { r -> record(id: r.id, meta: r.meta, gff: r.gunzip) }
        } else {
            ch_gff = channel.value(record(id: 'gff', meta: [:], gff: file(gff, checkIfExists: true)))
        }
        ch_gtf = GFFREAD(ch_gff, ch_no_path).map { r -> r.gtf }
    } else {
        ch_gtf = ch_no_path
    }

    //-------------------------------------
    // 2) Check if we actually have a FASTA
    //-------------------------------------
    if (fasta_provided && fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA(record(id: 'fasta', meta: [:], archive: file(fasta, checkIfExists: true))).map { r -> r.gunzip }
    } else if (fasta_provided) {
        ch_fasta = channel.value(file(fasta, checkIfExists: true))
    } else {
        ch_fasta = ch_no_path
    }

    //----------------------------------------
    // 3) Filter GTF if needed & FASTA present
    //----------------------------------------
    def filter_gtf_needed = (
        (!skip_alignment && aligner) ||
        (!skip_pseudo_alignment && pseudo_aligner) ||
        (!transcript_fasta)
    ) && !skip_gtf_filter

    // Set at the first step that supersedes ch_gtf (if any).
    ch_gtf_pre_filter = ch_no_path
    if (filter_gtf_needed && has_gtf) {
        ch_gtf_pre_filter = ch_gtf
        ch_gtf_filtered = CUSTOM_GTFFILTER(
            ch_gtf.map { item -> record(id: 'gtf', meta: [id: item.baseName + '.filtered'], gtf: item) },
            ch_fasta
        )
        ch_gtf = ch_gtf_filtered.map { r -> r.gtf }
    }

    //---------------------------------------------------
    // 4) Concatenate additional FASTA (if both are given)
    //---------------------------------------------------
    ch_additional_fasta_uncompressed = ch_no_path
    ch_fasta_pre_concat = ch_no_path
    ch_gtf_pre_concat   = ch_no_path
    if (fasta_provided && additional_fasta && has_gtf) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA(record(id: 'additional_fasta', meta: [:], archive: file(additional_fasta, checkIfExists: true))).map { r -> r.gunzip }
            ch_additional_fasta_uncompressed = ch_add_fasta
        } else {
            ch_add_fasta = channel.value(file(additional_fasta, checkIfExists: true))
        }

        ch_fasta_pre_concat = ch_fasta
        // Without a filter step the pre-concat GTF is the pre-filter GTF; record it once.
        if (filter_gtf_needed) {
            ch_gtf_pre_concat = ch_gtf
        } else {
            ch_gtf_pre_filter = ch_gtf
        }

        ch_catfasta = CUSTOM_CATADDITIONALFASTA(
            ch_fasta
                .combine(ch_gtf)
                .map { fasta_file, gtf_file -> record(id: 'genome_transcriptome', meta: [id: 'genome_transcriptome'], fasta: fasta_file, gtf: gtf_file) },
            ch_add_fasta,
            gencode ? "gene_type" : featurecounts_group_type
        )
        ch_fasta = ch_catfasta.map { r -> r.fasta }
        ch_gtf   = ch_catfasta.map { r -> r.gtf }
    } else if (fasta_provided && additional_fasta) {
        // No GTF to concatenate against: the additional FASTA is only uncompressed for publishing.
        ch_fasta_pre_concat = ch_fasta
        if (additional_fasta.endsWith('.gz')) {
            ch_additional_fasta_uncompressed = GUNZIP_ADDITIONAL_FASTA(record(id: 'additional_fasta', meta: [:], archive: file(additional_fasta, checkIfExists: true))).map { r -> r.gunzip }
        }
    }

    //------------------------------------------------------
    // 5) Uncompress gene BED or create from GTF if not given
    //------------------------------------------------------
    if (gene_bed && gene_bed.endsWith('.gz')) {
        ch_gene_bed = GUNZIP_GENE_BED(record(id: 'gene_bed', meta: [:], archive: file(gene_bed, checkIfExists: true))).map { r -> r.gunzip }
    } else if (gene_bed) {
        ch_gene_bed = channel.value(file(gene_bed, checkIfExists: true))
    } else if (prokaryotic && has_gtf) {
        // Prokaryotic annotations describe genes as CDS features, not exons, so
        // ea-utils/gtf2bed (which only reads `exon` rows) emits an empty BED.
        // gffread --bed derives intervals from any feature type.
        ch_gene_bed = GFFREAD_GENE_BED(
            ch_gtf.map { item -> record(id: item.baseName, meta: [id: item.baseName], gff: item) },
            ch_no_path
        ).map { r -> r.bed }
    } else if (has_gtf) {
        ch_gene_bed = EAUTILS_GTF2BED(ch_gtf.map { item -> record(id: item.baseName, meta: [id: item.baseName], gtf: item) }).map { r -> r.bed }
    } else {
        ch_gene_bed = ch_no_path
    }

    //----------------------------------------------------------------------
    // 6) Transcript FASTA:
    //    - If provided, decompress (optionally preprocess if GENCODE)
    //    - If not provided but have genome+GTF, create from them
    //----------------------------------------------------------------------
    ch_transcript_fasta_pre_gencode = ch_no_path
    ch_transcript_fasta_rsem_dir = ch_no_path
    if (transcript_fasta) {
        // Use user-provided transcript FASTA
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta_supplied = GUNZIP_TRANSCRIPT_FASTA(record(id: 'transcript_fasta', meta: [:], archive: file(transcript_fasta, checkIfExists: true))).map { r -> r.gunzip }
        } else {
            ch_transcript_fasta_supplied = channel.value(file(transcript_fasta, checkIfExists: true))
        }
        if (gencode) {
            ch_transcript_fasta_pre_gencode = ch_transcript_fasta_supplied
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(
                ch_transcript_fasta_supplied.map { fasta_file -> record(id: 'transcript_fasta', meta: [:], fasta: fasta_file) }
            ).map { r -> r.fasta }
        } else {
            ch_transcript_fasta = ch_transcript_fasta_supplied
        }
    } else if (fasta_provided && has_gtf && gffread_transcript_fasta) {
        // Use gffread to extract transcripts instead of RSEM
        // gffread handles CDS features correctly (e.g., prokaryotic annotations lack exon features)
        ch_transcript_fasta = GFFREAD_TRANSCRIPTS(
            ch_gtf.map { gtf_file -> record(id: 'transcripts', meta: [id: 'transcripts'], gff: gtf_file) },
            ch_fasta
        ).map { r -> r.gffread_fasta }
    } else if (fasta_provided && has_gtf && use_sentieon_star) {
        // Build transcripts from genome if we have it
        ch_rsem_reference = SENTIEON_MAKE_TRANSCRIPTS_FASTA(
            ch_fasta.map { fasta_file -> record(id: 'genome', meta: [id: 'genome'], fasta: fasta_file) }.combine(gtf: ch_gtf)
        )
        ch_transcript_fasta          = ch_rsem_reference.map { r -> r.transcript_fasta }
        ch_transcript_fasta_rsem_dir = ch_rsem_reference.map { r -> r.index } // unused here; published via the genome record's transcript_fasta_rsem_dir field
    } else if (fasta_provided && has_gtf) {
        // Build transcripts from genome if we have it
        ch_rsem_reference = MAKE_TRANSCRIPTS_FASTA(
            ch_fasta.map { fasta_file -> record(id: 'genome', meta: [id: 'genome'], fasta: fasta_file) }.combine(gtf: ch_gtf)
        )
        ch_transcript_fasta          = ch_rsem_reference.map { r -> r.transcript_fasta }
        ch_transcript_fasta_rsem_dir = ch_rsem_reference.map { r -> r.index } // unused here; published via the genome record's transcript_fasta_rsem_dir field
    } else {
        ch_transcript_fasta = ch_no_path
    }

    //-------------------------------------------------------
    // 7) FAI / chrom.sizes only if we actually have a genome
    //-------------------------------------------------------
    if (fasta_provided) {
        ch_faidx       = SAMTOOLS_FAIDX(ch_fasta.map { item -> record(id: 'genome', meta: [:], fasta: item, fai: null) }, true)
        ch_chrom_sizes = ch_faidx.map { r -> r.sizes }
        ch_fai         = ch_faidx.map { r -> r.fai }
    } else {
        ch_chrom_sizes = ch_no_path
        ch_fai         = ch_no_path
    }

    //-------------------------------------------------------------
    // 8) rRNA fastas (used by sortmerna index and bowtie2 rRNA removal)
    //-------------------------------------------------------------
    ch_rrna_fastas = channel.empty()

    // Load rRNA FASTAs when using sortmerna or bowtie2 for rRNA removal.
    // SortMeRNA's --ref option rejects gzipped FASTAs, so any .gz entries in the
    // manifest are decompressed first (the SortMeRNA v4.3 databases ship as .fasta.gz).
    if (ribo_removal_tool in ['sortmerna', 'bowtie2']) {
        def ribo_db = file(sortmerna_fasta_list)
        ch_rrna_inputs = channel.fromList(ribo_db.readLines())
            .map { row -> file(row) }

        ch_rrna_gunzipped = GUNZIP_RRNA_FASTAS(
            ch_rrna_inputs
                .filter { rrna_fasta -> rrna_fasta.name.endsWith('.gz') }
                .map { rrna_fasta -> record(id: rrna_fasta.name, meta: [:], archive: rrna_fasta) }
        )

        ch_rrna_fastas = ch_rrna_gunzipped
            .map { r -> r.gunzip }
            .mix(ch_rrna_inputs.filter { rrna_fasta -> !rrna_fasta.name.endsWith('.gz') })
    }

    //---------------------------------------------------------
    // 9) Kraken2 database (for contaminant screening)
    //---------------------------------------------------------
    if (contaminant_screening && kraken_db && kraken_db.endsWith('.tar.gz')) {
        ch_kraken_db = UNTAR_KRAKEN_DB(record(id: 'kraken_db', meta: [:], archive: file(kraken_db, checkIfExists: true))).map { r -> r.untar }
    } else if (contaminant_screening && kraken_db) {
        ch_kraken_db = channel.value(file(kraken_db, checkIfExists: true))
    } else {
        ch_kraken_db = ch_no_path
    }

    //---------------------------------------------------------
    // 10) Streamed reference/intermediate artifacts, one record per producer
    //---------------------------------------------------------
    // Each field is an independent, optionally user-supplied producer with no shared
    // key, so they stream as tagged records via mix() rather than fusing into one row.
    ch_references = ch_fasta.flatMap { p -> [ record(kind: 'fasta', file: p) ] }
        .mix(ch_fai.flatMap              { p -> [ record(kind: 'fai', file: p) ] })
        .mix(ch_gtf.flatMap              { p -> [ record(kind: 'gtf', file: p) ] })
        .mix(ch_gene_bed.flatMap         { p -> [ record(kind: 'gene_bed', file: p) ] })
        .mix(ch_transcript_fasta.flatMap { p -> [ record(kind: 'transcript_fasta', file: p) ] })
        .mix(ch_chrom_sizes.flatMap      { p -> [ record(kind: 'chrom_sizes', file: p) ] })
        .mix(ch_rrna_fastas.map          { f -> record(kind: 'rrna_fasta', file: f) })
        .mix(ch_kraken_db.flatMap        { p -> [ record(kind: 'kraken_db', file: p) ] })
        .filter { r -> taskOutputOrNull(r.file) != null }

    ch_intermediates = ch_gff_uncompressed.flatMap        { f -> [ record(kind: 'gff', file: f) ] }
        .mix(ch_additional_fasta_uncompressed.flatMap     { f -> [ record(kind: 'additional_fasta', file: f) ] })
        .mix(ch_gtf_pre_filter.flatMap                    { f -> [ record(kind: 'gtf_pre_filter', file: f) ] })
        .mix(ch_fasta_pre_concat.flatMap                  { f -> [ record(kind: 'fasta_pre_concat', file: f) ] })
        .mix(ch_gtf_pre_concat.flatMap                    { f -> [ record(kind: 'gtf_pre_concat', file: f) ] })
        .mix(ch_transcript_fasta_pre_gencode.flatMap      { f -> [ record(kind: 'tx_pre_gencode', file: f) ] })
        .mix(ch_transcript_fasta_rsem_dir.flatMap         { f -> [ record(kind: 'tx_rsem_dir', file: f) ] })
        .filter { r -> taskOutputOrNull(r.file) != null }

    // The genome streams stay separate named emits: RNASEQ and the output block consume each
    // artifact independently, and no per-sample key exists to fuse them into one record.
    emit:
    fasta:            Value<Path?> = ch_fasta                            // genome.fasta, null when absent
    fai:              Value<Path?> = ch_fai                              // genome.fai, null when no FASTA is given
    gtf:              Value<Path?> = ch_gtf                              // genome.gtf, null when absent
    gene_bed:         Value<Path?> = ch_gene_bed                         // gene.bed, null when absent
    transcript_fasta: Value<Path?> = ch_transcript_fasta                 // transcript.fasta, null when absent
    chrom_sizes:      Value<Path?> = ch_chrom_sizes                      // genome.sizes, null when absent
    rrna_fastas:      Channel<Path> = ch_rrna_fastas                     // rRNA fastas
    kraken_db:        Value<Path?> = ch_kraken_db                        // kraken2/db/, null when absent
    references:       Channel<GenomeArtifact> = ch_references            // one record per top-level reference file actually built or supplied
    intermediates:    Channel<GenomeArtifact> = ch_intermediates         // one record per superseded/incidental reference file
}
