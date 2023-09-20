//
// Prepare reference genome files
//

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { GFFREAD                           } from '../../modules/nf-core/gffread/main'
include { BBMAP_BBSPLIT                     } from '../../modules/nf-core/bbmap/bbsplit/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'
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
    ch_fasta              //      file: /path/to/genome.fasta
    ch_gtf                //      file: /path/to/genome.gtf
    ch_gff                //      file: /path/to/genome.gff
    ch_additional_fasta   //      file: /path/to/additional.fasta
    ch_transcript_fasta   //      file: /path/to/transcript.fasta
    ch_gene_bed           //      file: /path/to/gene.bed
    ch_splicesites        //      file: /path/to/splicesites.txt
    bbsplit_fasta_list    //      file: /path/to/bbsplit_fasta_list.txt
    star_index            // directory: /path/to/star/index/
    rsem_index            // directory: /path/to/rsem/index/
    salmon_index          // directory: /path/to/salmon/index/
    ch_hisat2_index       // directory: /path/to/hisat2/index/
    bbsplit_index         // directory: /path/to/rsem/index/
    gencode               //   boolean: whether the genome is from GENCODE
    is_aws_igenome        //   boolean: whether the genome files are from AWS iGenomes
    biotype               //    string: if additional fasta file is provided biotype value to use when appending entries to GTF file
    prepare_tool_indices  //      list: tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    //
    // Create GTF annotation from GFF3 if required
    //
    if (!ch_gtf && ch_gff) {
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Concatenate additional fasta file with reference fasta and gtf files
    //
    if (ch_additional_fasta) {
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_additional_fasta, biotype )

        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }

    //
    // Create transcript fasta file if required
    //
    if (ch_transcript_fasta) {
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
    // Create gene BED annotation file from GTF if required
    //
    if (!ch_gene_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Generate BBSplit index from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices && !bbsplit_index) {
        ch_bbsplit_fasta_list = Channel
            .from(file(bbsplit_fasta_list))
            .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
            .flatMap { id, fasta -> [ [ 'id', id ], [ 'fasta', file(fasta, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
            .groupTuple()
            .map { it -> it[1] } // Get rid of keys and keep grouped values
            .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module

        ch_bbsplit_index = BBMAP_BBSPLIT ( [ [:], [] ], [], ch_fasta, ch_bbsplit_fasta_list, true ).index
        ch_versions      = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star_salmon' in prepare_tool_indices && !star_index) {
        if (is_aws_igenome) {
            ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
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
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
        ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (!salmon_index && 'salmon' in prepare_tool_indices) {
        ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
        ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
