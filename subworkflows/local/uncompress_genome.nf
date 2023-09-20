//
// Uncompress reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/gunzip/main'

include { HISAT2_EXTRACTSPLICESITES         } from '../../modules/nf-core/hisat2/extractsplicesites/main'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/untar/main'

workflow UNCOMPRESS_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    gff                  //      file: /path/to/genome.gff
    additional_fasta     //      file: /path/to/additional.fasta
    transcript_fasta     //      file: /path/to/transcript.fasta
    gene_bed             //      file: /path/to/gene.bed
    splicesites          //      file: /path/to/splicesites.txt
    star_index           // directory: /path/to/star/index/
    rsem_index           // directory: /path/to/rsem/index/
    salmon_index         // directory: /path/to/salmon/index/
    hisat2_index         // directory: /path/to/hisat2/index/
    bbsplit_index        // directory: /path/to/rsem/index/
    prepare_tool_indices //      list: tools to uncompress indices for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }

    //
    // Uncompress GTF annotation file
    //
    ch_gtf = Channel.empty()
    if (gtf) {
        ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    }

    //
    // Uncompress GFF annotation file
    //
    ch_gff = Channel.empty()
    if (gff) {
        ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
    }

    //
    // Uncompress additional fasta file
    //
    ch_additional_fasta = Channel.empty()
    if (additional_fasta) {
        ch_additional_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], additional_fasta ] ).gunzip.map { it[1] }
        ch_versions         = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
    }

    //
    // Uncompress transcript fasta file
    //
    ch_transcript_fasta = Channel.empty()
    if (transcript_fasta) {
        ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
        ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file
    //
    ch_gene_bed = Channel.empty()
    if (gene_bed) {
        ch_gene_bed = GUNZIP_GENE_BED ( [ [:], gene_bed ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
    }

    //
    // Extract Splice Sites
    //
    ch_splicesites = Channel.empty()
    if ('hisat2' in prepare_tool_indices && !splicesites) {
        ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
        ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    }

    //
    // Uncompress STAR index
    //
    ch_star_index = Channel.empty()
    if ('star_salmon' in prepare_tool_indices && star_index) {
        ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
        ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
    }

    //
    // Uncompress RSEM index
    //
    ch_rsem_index = Channel.empty()
    if ('star_rsem' in prepare_tool_indices && rsem_index) {
        ch_rsem_index = UNTAR_RSEM_INDEX ( [ [:], rsem_index ] ).untar.map { it[1] }
        ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
    }

    //
    // Uncompress Salmon index
    //
    ch_salmon_index = Channel.empty()
    if ('salmon' in prepare_tool_indices && salmon_index) {
        ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
        ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
    }

    //
    // Uncompress HISAT2 index
    //
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices && hisat2_index) {
        ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
        ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
    }

    //
    // Uncompress BBSplit index
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices && bbsplit_index) {
        ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ [:], bbsplit_index ] ).untar.map { it[1] }
        ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    gff              = ch_gff                    // channel: path(genome.gff)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    additional_fasta = ch_additional_fasta       // channel: path(additional.fasta)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
