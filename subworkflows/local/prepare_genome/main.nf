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

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
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
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'

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
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  // boolean: whether the genome is from GENCODE
    featurecounts_group_type // string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  // string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'
    pseudo_aligner           // string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_gtf_filter          // boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             // boolean: Skip BBSplit for removal of non-reference genome reads
    skip_sortmerna           // boolean: Skip sortmerna for removal of reads mapping to sequences in sortmerna_fasta_list
    skip_alignment           // boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    // boolean: Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star             // boolean: whether to use sentieon STAR version

    main:
    // Versions collector
    ch_versions = Channel.empty()

    //---------------------------
    // 1) Uncompress GTF or GFF -> GTF
    //---------------------------
    ch_gtf = Channel.empty()
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ([ [:], file(gtf, checkIfExists: true) ]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        def ch_gff
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ([ [:], file(gff, checkIfExists: true) ]).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
        }
        ch_gtf      = GFFREAD(ch_gff, []).gtf.map { it[1] }
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //-------------------------------------
    // 2) Check if we actually have a FASTA
    //-------------------------------------
    def fasta_provided = (fasta ? true : false)

    ch_fasta = Channel.of([])
    if (fasta_provided) {
        // Uncompress FASTA if needed
        if (fasta.endsWith('.gz')) {
            ch_fasta    = GUNZIP_FASTA ([ [:], file(fasta, checkIfExists: true) ]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = Channel.value(file(fasta, checkIfExists: true))
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
        ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
    }

    //---------------------------------------------------
    // 4) Concatenate additional FASTA (if both are given)
    //---------------------------------------------------
    ch_add_fasta = Channel.empty()
    if (fasta_provided && additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA([ [:], file(additional_fasta, checkIfExists: true) ]).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta, checkIfExists: true))
        }

        CUSTOM_CATADDITIONALFASTA(
            ch_fasta.combine(ch_gtf).map { fasta, gtf -> [ [:], fasta, gtf ] },
            ch_add_fasta.map { [ [:], it ] },
            gencode ? "gene_type" : featurecounts_group_type
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map { it[1] }.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map { it[1] }.first()
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    //------------------------------------------------------
    // 5) Uncompress gene BED or create from GTF if not given
    //------------------------------------------------------
    ch_gene_bed = Channel.empty()
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ([ [:], file(gene_bed, checkIfExists: true) ]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed, checkIfExists: true))
        }
    } else {
        ch_gene_bed = GTF2BED(ch_gtf).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //----------------------------------------------------------------------
    // 6) Transcript FASTA:
    //    - If provided, decompress (optionally preprocess if GENCODE)
    //    - If not provided but have genome+GTF, create from them
    //----------------------------------------------------------------------
    ch_transcript_fasta = Channel.empty()
    if (transcript_fasta) {
        // Use user-provided transcript FASTA
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ([ [:], file(transcript_fasta, checkIfExists: true) ]).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta, checkIfExists: true))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else if (fasta_provided) {

        if(use_sentieon_star){
            // Build transcripts from genome if we have it
            ch_transcript_fasta = SENTIEON_MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta
            ch_versions         = ch_versions.mix(SENTIEON_MAKE_TRANSCRIPTS_FASTA.out.versions)
        } else {
            // Build transcripts from genome if we have it
            ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta
            ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
        }

    }

    //-------------------------------------------------------
    // 7) FAI / chrom.sizes only if we actually have a genome
    //-------------------------------------------------------
    ch_fai         = Channel.empty()
    ch_chrom_sizes = Channel.empty()
    if (fasta_provided) {
        CUSTOM_GETCHROMSIZES(ch_fasta.map { [ [:], it ] })
        ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    }

    //------------------------------------------------
    // 8) Determine which indices we actually want built
    //------------------------------------------------
    def prepare_tool_indices = []
    if (!skip_bbsplit)                                           { prepare_tool_indices << 'bbsplit' }
    if (!skip_sortmerna)                                         { prepare_tool_indices << 'sortmerna' }
    if ((!skip_alignment && aligner) || aligner == 'star_rsem')  { prepare_tool_indices << aligner }
    if (!skip_pseudo_alignment && pseudo_aligner)                { prepare_tool_indices << pseudo_aligner }

    //---------------------------------------------------------
    // 9) BBSplit index: uses FASTA only if we generate from scratch
    //---------------------------------------------------------
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            // Use user-provided bbsplit index
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ([ [:], file(bbsplit_index, checkIfExists: true) ]).untar.map { it[1] }
                ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            // Build it from scratch if we have FASTA
            Channel
                .from(file(bbsplit_fasta_list, checkIfExists: true))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fafile -> [ [ 'id', id ], [ 'fasta', file(fafile, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { it -> it[1] } // Get rid of keys and keep grouped values
                .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT(
                [ [:], [] ],
                [],
                ch_fasta,
                ch_bbsplit_fasta_list,
                true
            ).index
            ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
        }
        // else: no FASTA and no user-provided index -> remains empty
    }

    //-------------------------------------------------------------
    // 10) SortMeRNA index does not require the genome FASTA at all
    //-------------------------------------------------------------
    ch_sortmerna_index = Channel.empty()
    ch_rrna_fastas     = Channel.empty()
    if ('sortmerna' in prepare_tool_indices) {
        // We always need the rRNA FASTAs
        def ribo_db = file(sortmerna_fasta_list)
        ch_rrna_fastas = Channel.from(ribo_db.readLines())
            .map { row -> file(row) }

        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ([ [:], file(sortmerna_index, checkIfExists: true) ]).untar.map { it[1] }
                ch_versions        = ch_versions.mix(UNTAR_SORTMERNA_INDEX.out.versions)
            } else {
                ch_sortmerna_index = Channel.value([ [:], file(sortmerna_index, checkIfExists: true) ])
            }
        } else {
            // Build new SortMeRNA index from the rRNA references
            SORTMERNA_INDEX(
                Channel.of([ [], [] ]),
                ch_rrna_fastas.collect().map { [ 'rrna_refs', it ] },
                Channel.of([ [], [] ])
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
            ch_versions        = ch_versions.mix(SORTMERNA_INDEX.out.versions)
        }
    }

    //----------------------------------------------------
    // 11) STAR index (e.g. for 'star_salmon') -> needs FASTA if built
    //----------------------------------------------------
    ch_star_index = Channel.empty()
    if (prepare_tool_indices.intersect(['star_salmon', 'star_rsem'])) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ([ [:], file(star_index, checkIfExists: true) ]).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            // Build new STAR index
            // Possibly check AWS iGenome conditions
            def is_aws_igenome = false
            if (file(fasta, checkIfExists: true).getName() - '.gz' == 'genome.fa' && file(gtf, checkIfExists: true).getName() - '.gz' == 'genes.gtf') {
                is_aws_igenome = true
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES(ch_fasta, ch_gtf).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE(
                    ch_fasta.map { [ [:], it ] },
                    ch_gtf.map   { [ [:], it ] }
                ).index.map { it[1] }
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //------------------------------------------------
    // 12) RSEM index -> needs FASTA & GTF if built
    //------------------------------------------------
    ch_rsem_index = Channel.empty()
    if ('star_rsem' in prepare_tool_indices) {
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ([ [:], file(rsem_index, checkIfExists: true) ]).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
            } else {
                ch_rsem_index = Channel.value(file(rsem_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {

            if(use_sentieon_star){
                ch_rsem_index = SENTIEON_RSEM_PREPAREREFERENCE_GENOME(ch_fasta, ch_gtf).index
                ch_versions   = ch_versions.mix(SENTIEON_RSEM_PREPAREREFERENCE_GENOME.out.versions)
            }else{
                ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME(ch_fasta, ch_gtf).index
                ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
            }

        }
    }

    //---------------------------------------------------------
    // 13) HISAT2 index -> needs FASTA & GTF if built
    //---------------------------------------------------------
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        // splicesites
        if (splicesites) {
            ch_splicesites = Channel.value(file(splicesites, checkIfExists: true))
        }
        else if (fasta_provided) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES(ch_gtf.map { [ [:], it ] }).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        }
        // the index
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ([ [:], file(hisat2_index, checkIfExists: true) ]).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(hisat2_index, checkIfExists: true))
            }
        }
        else if (fasta_provided) {
            ch_hisat2_index = HISAT2_BUILD(
                ch_fasta.map { [ [:], it ] },
                ch_gtf.map   { [ [:], it ] },
                ch_splicesites.map { [ [:], it ] }
            ).index.map { it[1] }
            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    //------------------------------------------------------
    // 14) Salmon index -> can skip genome if transcript_fasta is enough
    //------------------------------------------------------

    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else if ('salmon' in prepare_tool_indices) {
        if (ch_transcript_fasta && fasta_provided) {
            // build from transcript FASTA + genome FASTA
            ch_salmon_index = SALMON_INDEX(ch_fasta, ch_transcript_fasta).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
        else if (ch_transcript_fasta) {
            // some Salmon module can run with just a transcript FASTA
            ch_salmon_index = SALMON_INDEX([], ch_transcript_fasta).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    //--------------------------------------------------
    // 15) Kallisto index -> only needs transcript FASTA
    //--------------------------------------------------
    ch_kallisto_index = Channel.empty()
    if (kallisto_index) {
        if (kallisto_index.endsWith('.tar.gz')) {
            ch_kallisto_index = UNTAR_KALLISTO_INDEX ( [ [:], kallisto_index ] ).untar
            ch_versions     = ch_versions.mix(UNTAR_KALLISTO_INDEX.out.versions)
        } else {
            ch_kallisto_index = Channel.value([[:], file(kallisto_index)])
        }
    } else {
        if ('kallisto' in prepare_tool_indices) {
            ch_kallisto_index = KALLISTO_INDEX ( ch_transcript_fasta.map { [ [:], it] } ).index
            ch_versions     = ch_versions.mix(KALLISTO_INDEX.out.versions)
        }
    }

    //------------------
    // 16) Emit channels
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
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
