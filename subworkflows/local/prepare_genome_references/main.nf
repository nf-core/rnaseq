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

workflow PREPARE_GENOME_REFERENCES {

    take:
    fasta                    // file: /path/to/genome.fasta (optional!)
    gtf                      // file: /path/to/genome.gtf
    gff                      // file: /path/to/genome.gff
    additional_fasta         // file: /path/to/additional.fasta
    transcript_fasta         // file: /path/to/transcript.fasta
    gene_bed                 // file: /path/to/gene.bed
    sortmerna_fasta_list     // file: /path/to/sortmerna_fasta_list.txt
    kraken_db                // path: /path/to/kraken2/db/ or .tar.gz archive
    gencode                  // boolean: whether the genome is from GENCODE
    gffread_transcript_fasta // boolean: use gffread instead of RSEM for transcript FASTA extraction
    featurecounts_group_type // string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  // string: Specifies the alignment algorithm to use
    pseudo_aligner           // string: Specifies the pseudo aligner to use
    skip_gtf_filter          // boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    ribo_removal_tool        // string: Tool for rRNA removal - 'sortmerna', 'ribodetector', or 'bowtie2' (null if skip)
    skip_alignment           // boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    // boolean: Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star        // boolean: whether to use sentieon STAR version
    contaminant_screening    // string: contaminant screening tool ('kraken2', 'kraken2_bracken', 'sylph', or null)
    prokaryotic              // boolean: whether the genome is prokaryotic (CDS-only annotation - use gffread --bed for gene BED since ea-utils/gtf2bed only handles exon features)

    main:
    //---------------------------
    // 1) Uncompress GTF or GFF -> GTF
    //---------------------------
    ch_gtf = channel.empty()
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ([ [:], file(gtf, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_gtf = channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        def ch_gff
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ([ [:], file(gff, checkIfExists: true) ]).gunzip
        } else {
            ch_gff = channel.value(file(gff, checkIfExists: true)).map { item -> [ [:], item ] }
        }
        ch_gtf      = GFFREAD(ch_gff, []).gtf.map { tuple -> tuple[1] }
    }

    //-------------------------------------
    // 2) Check if we actually have a FASTA
    //-------------------------------------
    def fasta_provided = (fasta ? true : false)

    ch_fasta = channel.of([])
    if (fasta_provided) {
        // Uncompress FASTA if needed
        if (fasta.endsWith('.gz')) {
            ch_fasta    = GUNZIP_FASTA ([ [:], file(fasta, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_fasta = channel.value(file(fasta, checkIfExists: true))
        }
    }

    //----------------------------------------
    // 3) Filter GTF if needed & FASTA present
    //----------------------------------------
    def filter_gtf_needed = (
        (!skip_alignment && aligner) ||
        (!skip_pseudo_alignment && pseudo_aligner) ||
        (!transcript_fasta)
    ) && !skip_gtf_filter

    if (filter_gtf_needed) {
        CUSTOM_GTFFILTER(
            ch_gtf.map { item -> [ [id: item.baseName + '.filtered'], item ] },
            fasta_provided
                ? ch_fasta.map { item -> [ [id: 'genome'], item ] }
                : channel.value([ [id: 'no_fasta'], [] ])
        )
        ch_gtf      = CUSTOM_GTFFILTER.out.gtf.map { _meta, filtered_gtf -> filtered_gtf }.first()
    }

    //---------------------------------------------------
    // 4) Concatenate additional FASTA (if both are given)
    //---------------------------------------------------
    ch_add_fasta = channel.empty()
    if (fasta_provided && additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA([ [:], file(additional_fasta, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_add_fasta = channel.value(file(additional_fasta, checkIfExists: true))
        }

        CUSTOM_CATADDITIONALFASTA(
            ch_fasta.combine(ch_gtf).map { fasta_file, gtf_file -> [ [id: 'genome_transcriptome'], fasta_file, gtf_file ] },
            ch_add_fasta.map { item -> [ [id: 'genome_transcriptome'], item ] },
            gencode ? "gene_type" : featurecounts_group_type
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map { tuple -> tuple[1] }.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map { tuple -> tuple[1] }.first()
    }

    //------------------------------------------------------
    // 5) Uncompress gene BED or create from GTF if not given
    //------------------------------------------------------
    ch_gene_bed = channel.empty()
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ([ [:], file(gene_bed, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_gene_bed = channel.value(file(gene_bed, checkIfExists: true))
        }
    } else if (prokaryotic) {
        // Prokaryotic annotations describe genes as CDS features, not exons, so
        // ea-utils/gtf2bed (which only reads `exon` rows) emits an empty BED.
        // gffread --bed derives intervals from any feature type.
        ch_gene_bed = GFFREAD_GENE_BED(
            ch_gtf.map { item -> [ [id: item.baseName], item ] },
            []
        ).bed.map { _meta, bed -> bed }
    } else {
        ch_gene_bed = EAUTILS_GTF2BED(ch_gtf.map { item -> [ [id: item.baseName], item ] }).bed.map { _meta, bed -> bed }
    }

    //----------------------------------------------------------------------
    // 6) Transcript FASTA:
    //    - If provided, decompress (optionally preprocess if GENCODE)
    //    - If not provided but have genome+GTF, create from them
    //----------------------------------------------------------------------
    ch_transcript_fasta = channel.empty()
    if (transcript_fasta) {
        // Use user-provided transcript FASTA
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ([ [:], file(transcript_fasta, checkIfExists: true) ]).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_transcript_fasta = channel.value(file(transcript_fasta, checkIfExists: true))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
        }
    } else if (fasta_provided) {

        if (gffread_transcript_fasta) {
            // Use gffread to extract transcripts instead of RSEM
            // gffread handles CDS features correctly (e.g., prokaryotic annotations lack exon features)
            GFFREAD_TRANSCRIPTS(
                ch_gtf.map { gtf_file -> [ [id: 'transcripts'], gtf_file ] },
                ch_fasta
            )
            ch_transcript_fasta = GFFREAD_TRANSCRIPTS.out.gffread_fasta.map { _meta, fasta_file -> fasta_file }
        } else if (use_sentieon_star) {
            // Build transcripts from genome if we have it
            ch_transcript_fasta = SENTIEON_MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta
        } else {
            // Build transcripts from genome if we have it
            ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta
        }

    }

    //-------------------------------------------------------
    // 7) FAI / chrom.sizes only if we actually have a genome
    //-------------------------------------------------------
    ch_fasta_fai   = channel.empty()
    ch_chrom_sizes = channel.empty()
    if (fasta_provided) {
        SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, true)
        ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { tuple -> tuple[1] }
        ch_fasta_fai   = ch_fasta
            .combine(SAMTOOLS_FAIDX.out.fai.map { _meta, fai_file -> fai_file })
            .map { fasta_file, fai_file -> [ [:], fasta_file, fai_file ] }
            .first()
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
        def ch_rrna_inputs = channel.from(ribo_db.readLines())
            .map { row -> file(row) }
            .branch { rrna_fasta ->
                gz:    rrna_fasta.name.endsWith('.gz')
                plain: true
            }

        ch_rrna_fastas = GUNZIP_RRNA_FASTAS(ch_rrna_inputs.gz.map { rrna_fasta -> [ [:], rrna_fasta ] })
            .gunzip
            .map { tuple -> tuple[1] }
            .mix(ch_rrna_inputs.plain)
    }

    //---------------------------------------------------------
    // 9) Kraken2 database (for contaminant screening)
    //---------------------------------------------------------
    ch_kraken_db = channel.empty()
    if (contaminant_screening && kraken_db) {
        if (kraken_db.endsWith('.tar.gz')) {
            ch_kraken_db = UNTAR_KRAKEN_DB ( [ [:], file(kraken_db, checkIfExists: true) ] ).untar.map { tuple -> tuple[1] }
        } else {
            ch_kraken_db = channel.value(file(kraken_db, checkIfExists: true))
        }
    }

    emit:
    fasta_fai        = ch_fasta_fai              // channel: [ meta, path(genome.fasta), path(genome.fai) ]
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    rrna_fastas      = ch_rrna_fastas            // channel: path(rrna_fastas)
    kraken_db        = ch_kraken_db              // channel: path(kraken2/db/)
}
