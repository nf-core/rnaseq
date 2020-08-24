/*
 * Map reads, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BWA_MEM           } from '../software/bwa/mem/main'
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'

workflow MAP_BWA_MEM {
    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_index         //    path: /path/to/index
    ch_fasta         //    path: /path/to/genome.fasta
    bwa_mem_options  //     map: options for BWA MEM module
    samtools_options //     map: options for SAMTools modules

    main:
    BWA_MEM(ch_reads, ch_index, ch_fasta, bwa_mem_options)
    BAM_SORT_SAMTOOLS(BWA_MEM.out.bam, samtools_options)

    emit:
    bam = BAM_SORT_SAMTOOLS.out.bam                           // channel: [ val(meta), [ bam ] ]
    bai = BAM_SORT_SAMTOOLS.out.bai                           // channel: [ val(meta), [ bai ] ]
    stats = BAM_SORT_SAMTOOLS.out.stats                       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat                 // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats                 // channel: [ val(meta), [ idxstats ] ]
    bwa_version = BWA_MEM.out.version                         //    path: *.version.txt
    samtools_version = BAM_SORT_SAMTOOLS.out.samtools_version //    path: *.version.txt
}
