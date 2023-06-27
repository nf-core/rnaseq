//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/gunzip/main'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { GFFREAD                           } from '../../modules/nf-core/gffread/main'
include { BBMAP_BBSPLIT                     } from '../../modules/nf-core/bbmap/bbsplit/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_EXTRACTSPLICESITES         } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                      } from '../../modules/nf-core/hisat2/build/main'
include { SALMON_INDEX                      } from '../../modules/nf-core/salmon/index/main'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../modules/nf-core/rsem/preparereference/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference/main'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../modules/local/gtf2bed'
include { CAT_ADDITIONAL_FASTA                 } from '../../modules/local/cat_additional_fasta'
include { GTF_GENE_FILTER                      } from '../../modules/local/gtf_gene_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../modules/local/star_genomegenerate_igenomes'

workflow PREPARE_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    gff                  //      file: /path/to/genome.gff
    additional_fasta     //      file: /path/to/additional.fasta
    transcript_fasta     //      file: /path/to/transcript.fasta
    gene_bed             //      file: /path/to/gene.bed
    splicesites          //      file: /path/to/splicesites.txt
    bbsplit_fasta_list   //      file: /path/to/bbsplit_fasta_list.txt
    star_index           // directory: /path/to/star/index/
    rsem_index           // directory: /path/to/rsem/index/
    salmon_index         // directory: /path/to/salmon/index/
    hisat2_index         // directory: /path/to/hisat2/index/ 
    bbsplit_index        // directory: /path/to/rsem/index/
    gencode              //   boolean: whether the genome is from GENCODE
    is_aws_igenome       //   boolean: whether the genome files are from AWS iGenomes
    biotype              //    string: if additional fasta file is provided biotype value to use when appending entries to GTF file
    prepare_tool_indices //      list: tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(gff))
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta))
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, biotype )
        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed))
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (gencode) { 
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_filter_gtf = GTF_GENE_FILTER ( ch_fasta, ch_gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ [:], bbsplit_index ] ).untar.map { it[1] }
                ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index))
            }
        } else {
            Channel
                .from(file(bbsplit_fasta_list))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fasta -> [ [ 'id', id ], [ 'fasta', file(fasta, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { it -> it[1] } // Get rid of keys and keep grouped values
                .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT ( [ [:], [] ], [], ch_fasta, ch_bbsplit_fasta_list, true ).index
            ch_versions      = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star_salmon' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
        } else {
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index = Channel.empty()
    if ('star_rsem' in prepare_tool_indices) {
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( [ [:], rsem_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
            } else {
                ch_rsem_index = Channel.value(file(rsem_index))
            }
        } else {
            ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.value(file(splicesites))
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(hisat2_index))
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
