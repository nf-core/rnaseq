/*
 * Alignment with HISAT2
 */

include { UNTAR                     } from '../process/untar'
include { HISAT2_EXTRACTSPLICESITES } from '../process/hisat2_extractsplicesites'
include { HISAT2_BUILD              } from '../process/hisat2_build'
include { HISAT2_ALIGN              } from '../process/hisat2_align'

workflow ALIGN_HISAT2 {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index         //    file: /path/to/star/index/
    fasta         //    file: /path/to/genome.fasta
    gtf           //    file: /path/to/genome.gtf
    splicesites   //    file: /path/to/genome.splicesites.txt
    index_options //     map: options for hisat2_build module
    align_options //     map: options for hisat2_align module

    main:
    /*
     * Uncompress HISAT2 index or generate from scratch if required
    */
    if (!splicesites) {
        ch_splicesites = HISAT2_EXTRACTSPLICESITES ( gtf, index_options ).txt
    } else {
        ch_splicesites = file(splicesites)
    }
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index, index_options ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = HISAT2_BUILD ( fasta, gtf, ch_splicesites, index_options ).index
    }

    /*
     * Map reads with HISAT2
     */
    HISAT2_ALIGN ( reads, ch_index, ch_splicesites, align_options )

    emit:
    bam     = HISAT2_ALIGN.out.bam     // channel: [ val(meta), bam   ]
    summary = HISAT2_ALIGN.out.summary // channel: [ val(meta), log   ]
    fastq   = HISAT2_ALIGN.out.fastq   // channel: [ val(meta), fastq ]
    version = HISAT2_ALIGN.out.version //    path: *.version.txt
}
