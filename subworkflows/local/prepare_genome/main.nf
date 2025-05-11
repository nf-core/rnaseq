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

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    gene_bed                 //      file: /path/to/gene.bed
    splicesites              //      file: /path/to/splicesites.txt
    bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    rsem_index               // directory: /path/to/rsem/index/
    salmon_index             // directory: /path/to/salmon/index/
    kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'
    pseudo_aligner           //    string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             //   boolean: Skip BBSplit for removal of non-reference genome reads
    skip_sortmerna           //   boolean: Skip sortmerna for removal of reads mapping to sequences in sortmerna_fasta_list
    skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    //   boolean: Skip all of the pseudoalignment-based processes within the pipeline

    main:
    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( file(fasta) )
    } else {
        ch_fasta = file(fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf = GUNZIP_GTF ( file(gtf) )
            } else {
                ch_gtf = file(gtf)
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff = GUNZIP_GFF ( file(gff) )
            } else {
                ch_gff = file(gff)
            }
            ch_gtf = GFFREAD ( ch_gff ).gtf
        }

        // Determine whether to filter the GTF or not
        def filter_gtf =
            ((
                // Condition 1: Alignment is required and aligner is set
                !skip_alignment && aligner
            ) ||
            (
                // Condition 2: Pseudoalignment is required and pseudoaligner is set
                !skip_pseudo_alignment && pseudo_aligner
            ) ||
            (
                // Condition 3: Transcript FASTA file is not provided
                !transcript_fasta
            )) &&
            (
                // Condition 4: --skip_gtf_filter is not provided
                !skip_gtf_filter
            )
        if (filter_gtf) {
            ch_gtf = GTF_FILTER ( ch_fasta, ch_gtf )
        }
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : featurecounts_group_type
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( file(additional_fasta) )
        } else {
            ch_add_fasta = file(additional_fasta)
        }

        def out = CUSTOM_CATADDITIONALFASTA (
            ch_fasta,
            ch_gtf,
            ch_add_fasta,
            biotype
        )
        ch_fasta    = out.fasta
        ch_gtf      = out.gtf
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( file(gene_bed) )
        } else {
            ch_gene_bed = file(gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf )
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( file(transcript_fasta) )
        } else {
            ch_transcript_fasta = file(transcript_fasta)
        }
        if (gencode) {
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
    }

    //
    // Create chromosome sizes file
    //
    chromsizes = CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_fai         = chromsizes.fai
    ch_chrom_sizes = chromsizes.sizes

    //
    // Get list of indices that need to be created
    //
    def prepare_tool_indices = []
    if (!skip_bbsplit) { prepare_tool_indices << 'bbsplit' }
    if (!skip_sortmerna) { prepare_tool_indices << 'sortmerna' }
    if (!skip_alignment) { prepare_tool_indices << aligner }
    if (!skip_pseudo_alignment && pseudo_aligner) { prepare_tool_indices << pseudo_aligner }

    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = null
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [:], file(bbsplit_index) )
            } else {
                ch_bbsplit_index = file(bbsplit_index)
            }
        } else {
            other_refs = file(bbsplit_fasta_list)
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .collectMany { id, fasta_ ->
                    [ [ 'id', id ], [ 'fasta', file(fasta_, checkIfExists: true) ] ]
                } // Flatten entries to be able to groupTuple by a common key
                .groupBy { _key, value -> value } // Get rid of keys and keep grouped values

            ch_bbsplit_index = BBMAP_BBSPLIT ( [:], [], null, ch_fasta, other_refs.id, other_refs.fasta, true ).index
        }
    }

    //
    // Uncompress sortmerna index or generate from scratch if required
    //
    ch_sortmerna_index = null
    if ('sortmerna' in prepare_tool_indices) {
        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ( [:], file(sortmerna_index) )
            } else {
                ch_sortmerna_index = file(sortmerna_index)
            }
        } else {
            ch_sortmerna_fastas = file(sortmerna_fasta_list)
                .readLines()
                .collect { row -> file(row, checkIfExists: true) }

            ch_sortmerna_index = SORTMERNA_INDEX ( [:], [], ch_sortmerna_fastas, null ).index
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = null
    if ('star_salmon' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [:], file(star_index) )
            } else {
                ch_star_index = file(star_index)
            }
        } else {
            // Check if an AWS iGenome has been provided to use the appropriate version of STAR
            def is_aws_igenome = false
            if (fasta && gtf) {
                if ((file(fasta).getName() - '.gz' == 'genome.fa') && (file(gtf).getName() - '.gz' == 'genes.gtf')) {
                    is_aws_igenome = true
                }
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf )
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf )
            }
        }
    }

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index = null
    if ('star_rsem' in prepare_tool_indices) {
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( [:], file(rsem_index) )
            } else {
                ch_rsem_index = file(rsem_index)
            }
        } else {
            ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME ( ch_fasta, ch_gtf ).index
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = null
    ch_hisat2_index = null
    if ('hisat2' in prepare_tool_indices) {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf )
        } else {
            ch_splicesites = file(splicesites)
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [:], file(hisat2_index) )
            } else {
                ch_hisat2_index = file(hisat2_index)
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta, ch_gtf, ch_splicesites )
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = null
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [:], file(salmon_index) )
        } else {
            ch_salmon_index = file(salmon_index)
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta )
        }
    }

    //
    // Uncompress Kallisto index or generate from scratch if required
    //
    ch_kallisto_index = null
    if (kallisto_index) {
        if (kallisto_index.endsWith('.tar.gz')) {
            ch_kallisto_index = UNTAR_KALLISTO_INDEX ( [:], file(kallisto_index) )
        } else {
            ch_kallisto_index = file(kallisto_index)
        }
    } else {
        if ('kallisto' in prepare_tool_indices) {
            ch_kallisto_index = KALLISTO_INDEX ( ch_transcript_fasta )
        }
    }

    emit:
    fasta            = ch_fasta                  // path(genome.fasta)
    gtf              = ch_gtf                    // path(genome.gtf)
    fai              = ch_fai                    // path(genome.fai)
    gene_bed         = ch_gene_bed               // path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // path(genome.sizes)
    splicesites      = ch_splicesites            // path(genome.splicesites.txt)
    bbsplit_index    = ch_bbsplit_index          // path(bbsplit/index/)
    sortmerna_index  = ch_sortmerna_index        // path(sortmerna/index/)
    star_index       = ch_star_index             // path(star/index/)
    rsem_index       = ch_rsem_index             // path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // path(hisat2/index/)
    salmon_index     = ch_salmon_index           // path(salmon/index/)
    kallisto_index   = ch_kallisto_index         // path(kallisto/index/)
}
