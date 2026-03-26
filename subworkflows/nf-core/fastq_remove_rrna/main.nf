include { BOWTIE2_ALIGN                               } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_PE             } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD                               } from '../../../modules/nf-core/bowtie2/build/main'
include { RIBODETECTOR                                } from '../../../modules/nf-core/ribodetector/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_BOWTIE2    } from '../../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BOWTIE2      } from '../../../modules/nf-core/samtools/view/main'
include { SEQKIT_REPLACE                              } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_U2T         } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_STATS                                } from '../../../modules/nf-core/seqkit/stats/main'
include { SORTMERNA                                 } from '../../../modules/nf-core/sortmerna/main'
include { SORTMERNA as SORTMERNA_INDEX              } from '../../../modules/nf-core/sortmerna/main'

//
// Function that parses seqkit stats TSV output to extract the mean read length
// for use with RiboDetector's -l parameter
//
def getReadLengthFromSeqkitStats(stats_file) {
    def lines = stats_file.text.readLines()
    if (lines.size() < 2) {
        return 100 // Default fallback
    }

    def header = lines[0].split('\t')
    def avgLenIdx = header.findIndexOf { col -> col == 'avg_len' }
    if (avgLenIdx < 0) {
        return 100 // Default fallback if column not found
    }

    // Calculate mean avg_len across all files in the stats output
    def avgLens = lines[1..-1].collect { line -> line.split('\t')[avgLenIdx] as float }
    def meanAvgLen = avgLens.sum() / avgLens.size()

    return Math.round(meanAvgLen) as int
}

workflow FASTQ_REMOVE_RRNA {
    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_rrna_fastas       // channel: one or more fasta files containing rrna sequences
    ch_sortmerna_index   // channel: /path/to/sortmerna/index/ (optional)
    ch_bowtie2_index     // channel: /path/to/bowtie2/index/ (optional)
    ribo_removal_tool    // string (enum): 'sortmerna', 'ribodetector', or 'bowtie2'
    make_sortmerna_index // boolean: Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index   // boolean: Whether to create a bowtie2 index before running bowtie2

    main:

    ch_multiqc_files = channel.empty()
    ch_filtered_reads = ch_reads

    // Individual output channels for workflow outputs
    ch_sortmerna_log     = channel.empty()
    ch_ribodetector_log  = channel.empty()
    ch_seqkit_stats      = channel.empty()
    ch_bowtie2_log       = channel.empty()
    ch_bowtie2_index_out = channel.empty()
    ch_seqkit_prefixed   = channel.empty()
    ch_seqkit_converted  = channel.empty()

    if (ribo_removal_tool == 'sortmerna') {
        ch_sortmerna_fastas = ch_rrna_fastas
            .collect()
            .map { fastas -> [[id: 'rrna_refs'], fastas] }

        if (make_sortmerna_index) {
            SORTMERNA_INDEX(
                [[], []],
                ch_sortmerna_fastas,
                [[], []],
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
        }

        SORTMERNA(
            ch_filtered_reads,
            ch_sortmerna_fastas,
            ch_sortmerna_index,
        )

        ch_filtered_reads = SORTMERNA.out.reads
        ch_sortmerna_log = SORTMERNA.out.log
        ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log)
    }
    else if (ribo_removal_tool == 'ribodetector') {
        // Run seqkit stats to determine average read length
        SEQKIT_STATS(
            ch_filtered_reads
        )

        ch_seqkit_stats = SEQKIT_STATS.out.stats
        ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats)

        // Join stats with reads and calculate read length for RiboDetector
        ch_filtered_reads
            .join(SEQKIT_STATS.out.stats)
            .multiMap { meta, reads, stats ->
                def readLength = getReadLengthFromSeqkitStats(stats)
                reads: [meta, reads]
                length: readLength
            }
            .set { ch_reads_with_length }

        RIBODETECTOR(
            ch_reads_with_length.reads,
            ch_reads_with_length.length,
        )

        ch_filtered_reads = RIBODETECTOR.out.fastq
        ch_ribodetector_log = RIBODETECTOR.out.log
        ch_multiqc_files = ch_multiqc_files.mix(RIBODETECTOR.out.log)
        // Note: ribodetector versions collected via topic
    }
    else if (ribo_removal_tool == 'bowtie2') {
        if (make_bowtie2_index) {
            // Process each rRNA file to add unique prefixes and convert U to T
            // This prevents duplicate sequence IDs in SAM header when combining databases
            ch_rrna_fastas
                .map { fasta_file -> [[id: fasta_file.baseName], fasta_file] }
                .set { ch_rrna_with_meta }

            // Step 1: Add filename prefixes to sequence headers
            SEQKIT_REPLACE(
                ch_rrna_with_meta
            )
            ch_seqkit_prefixed = SEQKIT_REPLACE.out.fastx

            // Step 2: Convert U to T in sequences (RNA to DNA)
            SEQKIT_REPLACE.out.fastx
                .map { meta, fasta_file -> [[id: "${meta.id}_dna"], fasta_file] }
                .set { ch_prefixed_fastas }

            SEQKIT_REPLACE_U2T(
                ch_prefixed_fastas
            )
            ch_seqkit_converted = SEQKIT_REPLACE_U2T.out.fastx

            // Collect processed files (already prefixed and U->T converted)
            SEQKIT_REPLACE_U2T.out.fastx
                .map { _meta, fasta_file -> fasta_file }
                .collectFile(name: 'rrna_combined_dna.fasta', newLine: true)
                .map { fasta_file -> [[id: 'rrna_refs'], fasta_file] }
                .set { ch_combined_fasta }

            BOWTIE2_BUILD(
                ch_combined_fasta
            )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index.first()
            ch_bowtie2_index_out = BOWTIE2_BUILD.out.index
        }

        // Branch reads by single-end vs paired-end for different filtering strategies
        ch_filtered_reads
            .branch { meta, _reads ->
                single_end: meta.single_end
                paired_end: !meta.single_end
            }
            .set { ch_reads_for_bowtie2 }

        // For single-end reads: bowtie2's --un-gz works correctly
        // save_unaligned=true outputs unmapped reads directly
        BOWTIE2_ALIGN(
            ch_reads_for_bowtie2.single_end,
            ch_bowtie2_index,
            [[], []],  // No reference fasta needed
            true,      // save_unaligned - for single-end this works correctly
            false,     // sort_bam - not needed
        )

        ch_bowtie2_log = BOWTIE2_ALIGN.out.log
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

        // For paired-end reads: bowtie2's --un-conc-gz outputs pairs that didn't
        // align concordantly, which INCLUDES pairs where one mate aligned.
        // We need to filter via samtools to get pairs where BOTH mates are unmapped.
        BOWTIE2_ALIGN_PE(
            ch_reads_for_bowtie2.paired_end,
            ch_bowtie2_index,
            [[], []],  // No reference fasta needed for BAM output
            false,     // save_unaligned - we'll extract from BAM instead
            false,     // sort_bam - not needed
        )

        ch_bowtie2_log = ch_bowtie2_log.mix(BOWTIE2_ALIGN_PE.out.log)
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN_PE.out.log)

        // Filter BAM for read pairs where BOTH mates are unmapped (flag 12 = 4 + 8)
        // This removes any pair where at least one mate aligned to rRNA
        SAMTOOLS_VIEW_BOWTIE2(
            BOWTIE2_ALIGN_PE.out.bam.map { meta, bam_file -> [meta, bam_file, []] },
            [[], [], []],  // No reference fasta
            [],        // No qname file
            []         // No index format
        )
        // Note: samtools/view versions collected via topic

        // Convert filtered BAM back to paired FASTQ
        SAMTOOLS_FASTQ_BOWTIE2(
            SAMTOOLS_VIEW_BOWTIE2.out.bam,
            false  // not interleaved
        )

        // Combine single-end and paired-end results
        BOWTIE2_ALIGN.out.fastq
            .mix(SAMTOOLS_FASTQ_BOWTIE2.out.fastq)
            .set { ch_filtered_reads }
    }

    emit:
    reads            = ch_filtered_reads   // channel: [ val(meta), [ reads ] ]
    multiqc_files    = ch_multiqc_files    // channel: [ val(meta), [ log files ] ]
    sortmerna_log    = ch_sortmerna_log    // channel: [ val(meta), [ log ] ]
    ribodetector_log = ch_ribodetector_log // channel: [ val(meta), [ log ] ]
    seqkit_stats     = ch_seqkit_stats     // channel: [ val(meta), [ stats ] ]
    bowtie2_log      = ch_bowtie2_log      // channel: [ val(meta), [ log ] ]
    bowtie2_index    = ch_bowtie2_index_out // channel: [ val(meta), [ index ] ]
    seqkit_prefixed  = ch_seqkit_prefixed  // channel: [ val(meta), [ fasta ] ]
    seqkit_converted = ch_seqkit_converted // channel: [ val(meta), [ fasta ] ]
}
