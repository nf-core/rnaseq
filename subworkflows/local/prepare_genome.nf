//
// Uncompress and prepare reference genome files
//

params.genome_options       = [:]
params.index_options        = [:]
params.gffread_options      = [:]
params.star_index_options   = [:]
params.rsem_index_options   = [:]
params.hisat2_index_options = [:]
params.salmon_index_options = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_TRANSCRIPT_FASTA
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/modules/gunzip/main'                    addParams( options: params.genome_options       )
include { UNTAR as UNTAR_STAR_INDEX   } from '../../modules/nf-core/modules/untar/main'                     addParams( options: params.star_index_options   )
include { UNTAR as UNTAR_RSEM_INDEX   } from '../../modules/nf-core/modules/untar/main'                     addParams( options: params.index_options        )
include { UNTAR as UNTAR_HISAT2_INDEX } from '../../modules/nf-core/modules/untar/main'                     addParams( options: params.hisat2_index_options )
include { UNTAR as UNTAR_SALMON_INDEX } from '../../modules/nf-core/modules/untar/main'                     addParams( options: params.index_options        )
include { GFFREAD                     } from '../../modules/nf-core/modules/gffread/main'                   addParams( options: params.gffread_options      )
include { STAR_GENOMEGENERATE         } from '../../modules/nf-core/modules/star/genomegenerate/main'       addParams( options: params.star_index_options   )
include { HISAT2_EXTRACTSPLICESITES   } from '../../modules/nf-core/modules/hisat2/extractsplicesites/main' addParams( options: params.hisat2_index_options )
include { HISAT2_BUILD                } from '../../modules/nf-core/modules/hisat2/build/main'              addParams( options: params.hisat2_index_options )
include { SALMON_INDEX                } from '../../modules/nf-core/modules/salmon/index/main'              addParams( options: params.salmon_index_options )
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE             } from '../../modules/nf-core/modules/rsem/preparereference/main' addParams( options: params.rsem_index_options )
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_TRANSCRIPTS } from '../../modules/nf-core/modules/rsem/preparereference/main' addParams( options: params.genome_options     )

include { GTF2BED                     } from '../../modules/local/gtf2bed'              addParams( options: params.genome_options )
include { CAT_ADDITIONAL_FASTA        } from '../../modules/local/cat_additional_fasta' addParams( options: params.genome_options )
include { GTF_GENE_FILTER             } from '../../modules/local/gtf_gene_filter'      addParams( options: params.genome_options )
include { GET_CHROM_SIZES             } from '../../modules/local/get_chrom_sizes'      addParams( options: params.genome_options )

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list  : tools to prepare indices for
    biotype              // string: if additional fasta file is provided
                        //          biotype value to use when appending entries to GTF file

    main:

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gffread_version = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_gffread_version = GFFREAD.out.version
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( params.additional_fasta ).gunzip
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, biotype )
        ch_fasta = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf   = CAT_ADDITIONAL_FASTA.out.gtf
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( params.gene_bed ).gunzip
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf )
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    ch_rsem_version = Channel.empty()
    if (params.transcript_fasta) {
        if (params.transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( params.transcript_fasta ).gunzip
        } else {
            ch_transcript_fasta = file(params.transcript_fasta)
        }
    } else {
        ch_transcript_fasta = RSEM_PREPAREREFERENCE_TRANSCRIPTS ( ch_fasta, GTF_GENE_FILTER ( ch_fasta, ch_gtf ) ).transcript_fasta
        ch_rsem_version     = RSEM_PREPAREREFERENCE_TRANSCRIPTS.out.version
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = GET_CHROM_SIZES ( ch_fasta ).sizes

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index   = Channel.empty()
    ch_star_version = Channel.empty()
    if ('star_salmon' in prepare_tool_indices) {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( params.star_index ).untar
            } else {
                ch_star_index = file(params.star_index)
            }
        } else {
            ch_star_index   = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_star_version = STAR_GENOMEGENERATE.out.version
        }
    }

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index   = Channel.empty()
    if ('star_rsem' in prepare_tool_indices) {
        if (params.rsem_index) {
            if (params.rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( params.rsem_index ).untar
            } else {
                ch_rsem_index = file(params.rsem_index)
            }
        } else {
            ch_rsem_index   = RSEM_PREPAREREFERENCE ( ch_fasta, ch_gtf ).index
            ch_rsem_version = RSEM_PREPAREREFERENCE.out.version
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites    = Channel.empty()
    ch_hisat2_index   = Channel.empty()
    ch_hisat2_version = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        if (!params.splicesites) {
            ch_splicesites    = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
            ch_hisat2_version = HISAT2_EXTRACTSPLICESITES.out.version
        } else {
            ch_splicesites = file(params.splicesites)
        }
        if (params.hisat2_index) {
            if (params.hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( params.hisat2_index ).untar
            } else {
                ch_hisat2_index = file(params.hisat2_index)
            }
        } else {
            ch_hisat2_index   = HISAT2_BUILD ( ch_fasta, ch_gtf, ch_splicesites ).index
            ch_hisat2_version = HISAT2_BUILD.out.version
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index   = Channel.empty()
    ch_salmon_version = Channel.empty()
    if ('salmon' in prepare_tool_indices) {
        if (params.salmon_index) {
            if (params.salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX ( params.salmon_index ).untar
            } else {
                ch_salmon_index = file(params.salmon_index)
            }
        } else {
            ch_salmon_index   = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
            ch_salmon_version = SALMON_INDEX.out.version
        }
    }

    emit:
    fasta            = ch_fasta            // path: genome.fasta
    gtf              = ch_gtf              // path: genome.gtf
    gene_bed         = ch_gene_bed         // path: gene.bed
    transcript_fasta = ch_transcript_fasta // path: transcript.fasta
    chrom_sizes      = ch_chrom_sizes      // path: genome.sizes
    splicesites      = ch_splicesites      // path: genome.splicesites.txt
    star_index       = ch_star_index       // path: star/index/
    rsem_index       = ch_rsem_index       // path: rsem/index/
    hisat2_index     = ch_hisat2_index     // path: hisat2/index/
    salmon_index     = ch_salmon_index     // path: salmon/index/
    star_version     = ch_star_version     // path: *.version.txt
    rsem_version     = ch_rsem_version     // path: *.version.txt
    hisat2_version   = ch_hisat2_version   // path: *.version.txt
    salmon_version   = ch_salmon_version   // path: *.version.txt
    gffread_version  = ch_gffread_version  // path: *.version.txt
}
