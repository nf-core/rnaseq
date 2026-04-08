//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_BOWTIE2_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KRAKEN_DB          } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { GFFREAD as GFFREAD_TRANSCRIPTS    } from '../../../modules/nf-core/gffread'
include { BOWTIE2_BUILD                     } from '../../../modules/nf-core/bowtie2/build'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_GENOMEGENERATE as PARABRICKS_STARGENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { KALLISTO_INDEX                    } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'
include { SENTIEON_RSEMPREPAREREFERENCE as SENTIEON_RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/sentieon/rsempreparereference'
include { SENTIEON_RSEMPREPAREREFERENCE as SENTIEON_MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/sentieon/rsempreparereference'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'

workflow PREPARE_GENOME {

    take:
    fasta                    // file: /path/to/genome.fasta (optional!)
    gtf                      // file: /path/to/genome.gtf
    gff                      // file: /path/to/genome.gff
    additional_fasta         // file: /path/to/additional.fasta
    transcript_fasta         // file: /path/to/transcript.fasta
    gene_bed                 // file: /path/to/gene.bed
    splicesites              // file: /path/to/splicesites.txt
    bbsplit_fasta_list       // file: /path/to/bbsplit_fasta_list.txt
    sortmerna_fasta_list     // file: /path/to/sortmerna_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    rsem_index               // directory: /path/to/rsem/index/
    salmon_index             // directory: /path/to/salmon/index/
    kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bowtie2_index            // directory: /path/to/bowtie2/index/
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    kraken_db                // path: /path/to/kraken2/db/ or .tar.gz archive
    gencode                  // boolean: whether the genome is from GENCODE
    gffread_transcript_fasta // boolean: use gffread instead of RSEM for transcript FASTA extraction
    featurecounts_group_type // string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  // string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem', 'hisat2', and 'bowtie2_salmon'
    pseudo_aligner           // string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_gtf_filter          // boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             // boolean: Skip BBSplit for removal of non-reference genome reads
    ribo_removal_tool        // string: Tool for rRNA removal - 'sortmerna', 'ribodetector', or 'bowtie2' (null if skip)
    skip_alignment           // boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    // boolean: Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star        // boolean: whether to use sentieon STAR version
    use_parabricks_star      // boolean: whether to use parabricks STAR version
    contaminant_screening    // string: contaminant screening tool ('kraken2', 'kraken2_bracken', 'sylph', or null)

    main:
    // Versions collector
    ch_versions = channel.empty()

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
        GTF_FILTER(ch_fasta, ch_gtf)
        ch_gtf      = GTF_FILTER.out.genome_gtf.first()
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
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
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
    } else {
        ch_gene_bed = GTF2BED(ch_gtf).bed
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
            ch_transcript_fasta = GFFREAD_TRANSCRIPTS.out.gffread_fasta.map { meta, fasta_file -> fasta_file }
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
    ch_fai         = channel.empty()
    ch_chrom_sizes = channel.empty()
    if (fasta_provided) {
        SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, true)
        ch_fai         = SAMTOOLS_FAIDX.out.fai.map { tuple -> tuple[1] }
        ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { tuple -> tuple[1] }
    }

    //------------------------------------------------
    // 8) Determine which indices we actually want built
    //------------------------------------------------
    def prepare_tool_indices = []
    if (!skip_bbsplit)                                           { prepare_tool_indices << 'bbsplit' }
    if (ribo_removal_tool == 'sortmerna')                        { prepare_tool_indices << 'sortmerna' }
    if ((!skip_alignment && aligner) || aligner == 'star_rsem')  { prepare_tool_indices << aligner }
    if (!skip_pseudo_alignment && pseudo_aligner)                { prepare_tool_indices << pseudo_aligner }

    //---------------------------------------------------------
    // 9) BBSplit index: uses FASTA only if we generate from scratch
    //---------------------------------------------------------
    ch_bbsplit_index = channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            // Use user-provided bbsplit index
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ([ [:], file(bbsplit_index, checkIfExists: true) ]).untar.map { tuple -> tuple[1] }
            } else {
                ch_bbsplit_index = channel.value(file(bbsplit_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            // Build it from scratch if we have FASTA
            channel
                .from(file(bbsplit_fasta_list, checkIfExists: true))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fafile -> [ [ 'id', id ], [ 'fasta', file(fafile, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { entry -> entry[1] } // Get rid of keys and keep grouped values
                .collect { item -> [ item ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT(
                [ [:], [] ],
                [],
                ch_fasta,
                ch_bbsplit_fasta_list,
                true
            ).index
        }
        // else: no FASTA and no user-provided index -> remains empty
    }

    //-------------------------------------------------------------
    // 10) rRNA fastas and SortMeRNA index
    //-------------------------------------------------------------
    ch_sortmerna_index = channel.empty()
    ch_rrna_fastas     = channel.empty()

    // Load rRNA FASTAs when using sortmerna or bowtie2 for rRNA removal
    if (ribo_removal_tool in ['sortmerna', 'bowtie2']) {
        def ribo_db = file(sortmerna_fasta_list)
        ch_rrna_fastas = channel.from(ribo_db.readLines())
            .map { row -> file(row) }
    }

    // Build SortMeRNA index only when using sortmerna
    if ('sortmerna' in prepare_tool_indices) {
        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ([ [:], file(sortmerna_index, checkIfExists: true) ]).untar.map { tuple -> tuple[1] }
            } else {
                ch_sortmerna_index = channel.value([ [:], file(sortmerna_index, checkIfExists: true) ])
            }
        } else {
            // Build new SortMeRNA index from the rRNA references
            SORTMERNA_INDEX(
                channel.of([ [], [] ]),
                ch_rrna_fastas.collect().map { refs -> [ 'rrna_refs', refs ] },
                channel.of([ [], [] ])
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
        }
    }

    //----------------------------------------------------
    // 11) STAR index (e.g. for 'star_salmon') -> needs FASTA if built
    //----------------------------------------------------
    ch_star_index = channel.empty()
    if (prepare_tool_indices.intersect(['star_salmon', 'star_rsem'])) {
        if (use_parabricks_star && fasta_provided) {
            // Parabricks needs its own STAR index built with its bundled STAR version
            ch_star_index = PARABRICKS_STARGENOMEGENERATE(
                ch_fasta.map { item -> [ [:], item ] },
                ch_gtf.map   { item -> [ [:], item ] }
            ).index.map { tuple -> tuple[1] }
        } else if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ([ [:], file(star_index, checkIfExists: true) ]).untar.map { tuple -> tuple[1] }
            } else {
                ch_star_index = channel.value(file(star_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            // Build new STAR index with current STAR version.
            // No need for iGenomes STAR 2.6.1d here - that's only for pre-built index compatibility.
            ch_star_index = STAR_GENOMEGENERATE(
                ch_fasta.map { item -> [ [:], item ] },
                ch_gtf.map { item -> [ [:], item ] }
            ).index.map { tuple -> tuple[1] }
        }
    }

    //------------------------------------------------
    // 12) RSEM index -> needs FASTA & GTF if built
    //------------------------------------------------
    ch_rsem_index = channel.empty()
    if ('star_rsem' in prepare_tool_indices) {
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ([ [:], file(rsem_index, checkIfExists: true) ]).untar.map { tuple -> tuple[1] }
            } else {
                ch_rsem_index = channel.value(file(rsem_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {

            if(use_sentieon_star){
                ch_rsem_index = SENTIEON_RSEM_PREPAREREFERENCE_GENOME(ch_fasta, ch_gtf).index
            }else{
                ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME(ch_fasta, ch_gtf).index
            }

        }
    }

    //---------------------------------------------------------
    // 13) HISAT2 index -> needs FASTA & GTF if built
    //---------------------------------------------------------
    ch_splicesites  = channel.empty()
    ch_hisat2_index = channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        // splicesites
        if (splicesites) {
            ch_splicesites = channel.value(file(splicesites, checkIfExists: true))
        }
        else if (fasta_provided) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES(ch_gtf.map { item -> [ [:], item ] }).txt.map { tuple -> tuple[1] }
        }
        // the index
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ([ [:], file(hisat2_index, checkIfExists: true) ]).untar.map { tuple -> tuple[1] }
            } else {
                ch_hisat2_index = channel.value(file(hisat2_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            ch_hisat2_index = HISAT2_BUILD(
                ch_fasta.map { item -> [ [:], item ] },
                ch_gtf.map { item -> [ [:], item ] },
                ch_splicesites.map { item -> [ [:], item ] }
            ).index.map { tuple -> tuple[1] }
        }
    }

    //---------------------------------------------------------
    // 14) Bowtie2 index -> built from transcript FASTA for Salmon alignment mode
    //---------------------------------------------------------
    ch_bowtie2_index = channel.empty()
    if ('bowtie2_salmon' in prepare_tool_indices) {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ([ [:], file(bowtie2_index, checkIfExists: true) ]).untar.map { meta, index -> index }
            } else {
                ch_bowtie2_index = channel.value(file(bowtie2_index, checkIfExists: true))
            }
        }
        else if (ch_transcript_fasta) {
            // Build Bowtie2 index from transcript FASTA for alignment-based Salmon quantification
            BOWTIE2_BUILD(
                ch_transcript_fasta.map { fasta_file -> [ [id: 'transcripts'], fasta_file ] }
            )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index.map { meta, index -> index }
        }
    }

    //------------------------------------------------------
    // 15) Salmon index -> can skip genome if transcript_fasta is enough
    //------------------------------------------------------

    ch_salmon_index = channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { tuple -> tuple[1] }
        } else {
            ch_salmon_index = channel.value(file(salmon_index))
        }
    } else if ('salmon' in prepare_tool_indices) {
        if (ch_transcript_fasta && fasta_provided) {
            ch_salmon_index = SALMON_INDEX(ch_fasta, ch_transcript_fasta).index
        }
        else if (ch_transcript_fasta) {
            ch_salmon_index = SALMON_INDEX([], ch_transcript_fasta).index
        }
    }

    //--------------------------------------------------
    // 16) Kallisto index -> only needs transcript FASTA
    //--------------------------------------------------
    ch_kallisto_index = channel.empty()
    if (kallisto_index) {
        if (kallisto_index.endsWith('.tar.gz')) {
            ch_kallisto_index = UNTAR_KALLISTO_INDEX ( [ [:], kallisto_index ] ).untar
        } else {
            ch_kallisto_index = channel.value([[:], file(kallisto_index)])
        }
    } else {
        if ('kallisto' in prepare_tool_indices) {
            ch_kallisto_index = KALLISTO_INDEX ( ch_transcript_fasta.map { item -> [ [:], item ] } ).index
        }
    }

    //---------------------------------------------------------
    // Kraken2 database (for contaminant screening)
    //---------------------------------------------------------
    ch_kraken_db = channel.empty()
    if (contaminant_screening && kraken_db) {
        if (kraken_db.endsWith('.tar.gz')) {
            ch_kraken_db = UNTAR_KRAKEN_DB ( [ [:], file(kraken_db, checkIfExists: true) ] ).untar.map { tuple -> tuple[1] }
        } else {
            ch_kraken_db = channel.value(file(kraken_db, checkIfExists: true))
        }
    }

    //------------------
    // 17) Emit channels
    //------------------

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    rrna_fastas      = ch_rrna_fastas            // channel: path(sortmerna_fasta_list)
    sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    bowtie2_index    = ch_bowtie2_index          // channel: path(bowtie2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
    kraken_db        = ch_kraken_db              // channel: path(kraken2/db/)
    versions         = ch_versions               // channel: [ versions.yml ]
}
