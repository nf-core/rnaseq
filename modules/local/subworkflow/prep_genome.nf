/*
 * Uncompress and prepare reference genome files
*/

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_BED12
    GUNZIP as GUNZIP_ADDITIONAL_FASTA
    GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../process/gunzip'
include { GFFREAD                     } from '../process/gffread'
include { GTF2BED                     } from '../process/gtf2bed'
include { CAT_ADDITIONAL_FASTA        } from '../process/cat_additional_fasta'

workflow PREP_GENOME {
    take:
    fasta             // file: /path/to/genome.fasta
    gtf               // file: /path/to/genome.gtf
    gff               // file: /path/to/genome.gff
    bed12             // file: /path/to/genome.bed12
    additional_fasta  // file: /path/to/additional.fasta
    options           //  map: module options for additional arguments and publishing files

    main:
    /*
     * PREPROCESSING: Uncompress genome fasta file if required
     */
    if (fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( fasta, options ).gunzip
    } else {
        ch_fasta = file(fasta)
    }

    /*
     * PREPROCESSING: Uncompress GTF annotation file or create from GFF3 if required
     */
    gffread_version = Channel.empty()
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( gtf, options ).gunzip
        } else {
            ch_gtf = file(gtf)
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( gff, options ).gunzip
        } else {
            ch_gff = file(gff)
        }
        ch_gtf = GFFREAD ( ch_gff, options ).gtf
        gffread_version = GFFREAD.out.version
    }

    /*
     * PREPROCESSING: Uncompress additional fasta file and concatenate with reference fasta and gtf files
     */
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( additional_fasta, options ).gunzip
        } else {
            ch_add_fasta = file(additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, options )
        ch_fasta = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf   = CAT_ADDITIONAL_FASTA.out.gtf
    }

    /*
     * PREPROCESSING: Uncompress BED12 annotation file or create from GTF if required
     */
    if (bed12) {
        if (bed12.endsWith('.gz')) {
            ch_bed12 = GUNZIP_BED12 ( bed12, options ).gunzip
        } else {
            ch_bed12 = file(bed12)
        }
    } else {
        ch_bed12 = GTF2BED ( ch_gtf, options )
    }

    emit:
    fasta           = ch_fasta // path: genome.fasta
    gtf             = ch_gtf   // path: genome.gtf
    bed12           = ch_bed12 // path: genome.bed12
    gffread_version            // path: *.version.txt
}
