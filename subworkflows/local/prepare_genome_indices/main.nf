//
// Build or load aligner / pseudo-aligner / filtering indices (STAR, RSEM, HISAT2, Bowtie2, Salmon, Kallisto, BBSplit, SortMeRNA)
//

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_BOWTIE2_INDEX      } from '../../../modules/nf-core/untar'

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
include { SENTIEON_RSEMPREPAREREFERENCE as SENTIEON_RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/sentieon/rsempreparereference'

include { STAR_GENOMEPARAMS_UPGRADE         } from '../../../modules/local/star_genomeparams_upgrade'

workflow PREPARE_GENOME_INDICES {

    take:
    ch_fasta_fai             // channel: [ meta, path(genome.fasta), path(genome.fai) ] - emitted from PREPARE_GENOME_REFERENCES
    ch_gtf                   // channel: path(genome.gtf) - emitted from PREPARE_GENOME_REFERENCES
    ch_transcript_fasta      // channel: path(transcript.fasta) - emitted from PREPARE_GENOME_REFERENCES
    ch_rrna_fastas           // channel: path(rrna_fastas) - emitted from PREPARE_GENOME_REFERENCES
    fasta_provided           // boolean: whether a genome FASTA was provided
    splicesites              // file: /path/to/splicesites.txt
    bbsplit_fasta_list       // file: /path/to/bbsplit_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    rsem_index               // directory: /path/to/rsem/index/
    salmon_index             // directory: /path/to/salmon/index/
    kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bowtie2_index            // directory: /path/to/bowtie2/index/
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    aligner                  // string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem', 'hisat2', and 'bowtie2_salmon'
    pseudo_aligner           // string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_bbsplit             // boolean: Skip BBSplit for removal of non-reference genome reads
    ribo_removal_tool        // string: Tool for rRNA removal - 'sortmerna', 'ribodetector', or 'bowtie2' (null if skip)
    skip_alignment           // boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    // boolean: Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star        // boolean: whether to use sentieon STAR version
    use_parabricks_star      // boolean: whether to use parabricks STAR version
    star_index_legacy        // boolean: whether the supplied star_index was built with STAR 2.6.x and needs genomeParameters.txt upgraded to the 2.7.4a metadata schema

    main:
    ch_fasta = ch_fasta_fai.map { _meta, fasta_file, _fai -> fasta_file }

    //------------------------------------------------
    // 1) Determine which indices we actually want built
    //------------------------------------------------
    def prepare_tool_indices = []
    if (!skip_bbsplit)                                           { prepare_tool_indices << 'bbsplit' }
    if (ribo_removal_tool == 'sortmerna')                        { prepare_tool_indices << 'sortmerna' }
    if ((!skip_alignment && aligner) || aligner == 'star_rsem')  { prepare_tool_indices << aligner }
    if (!skip_pseudo_alignment && pseudo_aligner)                { prepare_tool_indices << pseudo_aligner }

    //---------------------------------------------------------
    // 2) BBSplit index: uses FASTA only if we generate from scratch
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
    // 3) SortMeRNA index
    //-------------------------------------------------------------
    ch_sortmerna_index = channel.empty()

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
    // 4) STAR index (e.g. for 'star_salmon') -> needs FASTA if built
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
            // Pre-built STAR index supplied by the user. When star_index_legacy is set
            // (genomes-map opt-in for indices built with STAR 2.6.x, e.g. AWS iGenomes),
            // route through STAR_GENOMEPARAMS_UPGRADE to rewrite `versionGenome 20201` and
            // add the genomeType / genomeTransformType / genomeTransformVCF fields that
            // STAR 2.7.4a+ requires. Modern indices skip the adapter entirely.
            def ch_star_raw = star_index.endsWith('.tar.gz')
                ? UNTAR_STAR_INDEX([ [:], file(star_index, checkIfExists: true) ]).untar
                : channel.value([ [:], file(star_index, checkIfExists: true) ])
            ch_star_index = star_index_legacy
                ? STAR_GENOMEPARAMS_UPGRADE(ch_star_raw).index.map { tuple -> tuple[1] }
                : ch_star_raw.map { tuple -> tuple[1] }
        }
        else if (fasta_provided) {
            ch_star_index = STAR_GENOMEGENERATE(
                ch_fasta.map { item -> [ [:], item ] },
                ch_gtf.map { item -> [ [:], item ] }
            ).index.map { tuple -> tuple[1] }
        }
    }

    //------------------------------------------------
    // 5) RSEM index -> needs FASTA & GTF if built
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
    // 6) HISAT2 index -> needs FASTA & GTF if built
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
    // 7) Bowtie2 index -> built from transcript FASTA for Salmon alignment mode
    //---------------------------------------------------------
    ch_bowtie2_index = channel.empty()
    if ('bowtie2_salmon' in prepare_tool_indices) {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ([ [:], file(bowtie2_index, checkIfExists: true) ]).untar.map { _meta, index -> index }
            } else {
                ch_bowtie2_index = channel.value(file(bowtie2_index, checkIfExists: true))
            }
        }
        else if (ch_transcript_fasta) {
            // Build Bowtie2 index from transcript FASTA for alignment-based Salmon quantification
            BOWTIE2_BUILD(
                ch_transcript_fasta.map { fasta_file -> [ [id: 'transcripts'], fasta_file ] }
            )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index.map { _meta, index -> index }
        }
    }

    //------------------------------------------------------
    // 8) Salmon index -> can skip genome if transcript_fasta is enough
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
    // 9) Kallisto index -> only needs transcript FASTA
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

    emit:
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    bowtie2_index    = ch_bowtie2_index          // channel: path(bowtie2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
}
