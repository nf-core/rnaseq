nextflow.enable.types = true

include { BOWTIE2_ALIGN                            } from '../../../modules/nf-core/bowtie2/align'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_PE        } from '../../../modules/nf-core/bowtie2/align'
include { BOWTIE2_BUILD                            } from '../../../modules/nf-core/bowtie2/build'
include { RIBODETECTOR                             } from '../../../modules/nf-core/ribodetector'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_BOWTIE2 } from '../../../modules/nf-core/samtools/fastq'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BOWTIE2   } from '../../../modules/nf-core/samtools/view'
include { SEQKIT_REPLACE                           } from '../../../modules/nf-core/seqkit/replace'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_U2T     } from '../../../modules/nf-core/seqkit/replace'
include { SEQKIT_STATS                             } from '../../../modules/nf-core/seqkit/stats'
include { SORTMERNA                                } from '../../../modules/nf-core/sortmerna'
include { SORTMERNA as SORTMERNA_INDEX             } from '../../../modules/nf-core/sortmerna'
include { ReadsInput; RrnaReferences; FastqRemoveRrna } from '../../../modules/nf-core/types'
include { Bowtie2AlignResult } from '../../../modules/nf-core/bowtie2/align/main'
include { Bowtie2BuildResult } from '../../../modules/nf-core/bowtie2/build/main'
include { RibodetectorResult } from '../../../modules/nf-core/ribodetector/main'
include { SamtoolsFastqResult } from '../../../modules/nf-core/samtools/fastq/main'
include { SamtoolsViewResult } from '../../../modules/nf-core/samtools/view/main'
include { SeqkitReplaceResult } from '../../../modules/nf-core/seqkit/replace/main'
include { SeqkitStatsResult } from '../../../modules/nf-core/seqkit/stats/main'
include { SortmernaResult } from '../../../modules/nf-core/sortmerna/main'

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

record RrnaArgs {
    use_gpu_ribodetector: Boolean?
}

workflow FASTQ_REMOVE_RRNA {
    take:
    ch_reads: Channel<ReadsInput>
    ch_rrna_fastas: Channel<Path> // one or more fasta files containing rrna sequences
    ch_sortmerna_index: Value<Path> // sortmerna index (optional)
    ch_bowtie2_index: Value<Path> // bowtie2 index (optional)
    ribo_removal_tool: String // 'sortmerna', 'ribodetector', or 'bowtie2'
    make_sortmerna_index: Boolean // Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index: Boolean // Whether to create a bowtie2 index before running bowtie2
    tool_args: RrnaArgs // whether RiboDetector runs on a GPU

    main:

    // Run-level references built here, emitted separately from the per-sample results
    val_refs = channel.value(
        record(sortmerna_index: null, bowtie2_index: null, seqkit_prefixed: null, seqkit_converted: null)
    )

    if (ribo_removal_tool == 'sortmerna') {
        ch_sortmerna_fastas = ch_rrna_fastas.collect().map { refs -> refs.toList() }

        ch_sortmerna_idx = ch_sortmerna_index
        if (make_sortmerna_index) {
            def ch_sortmerna_built: Value<SortmernaResult> = SORTMERNA_INDEX(
                record(id: 'rrna_refs', meta: [:], reads: []),
                ch_sortmerna_fastas,
                channel.value(null as Path),
            )
            ch_sortmerna_idx = ch_sortmerna_built.map { r -> r.index }
            val_refs = ch_sortmerna_built.map { r ->
                record(sortmerna_index: r.index, bowtie2_index: null, seqkit_prefixed: null, seqkit_converted: null)
            }
        }

        def ch_sortmerna: Channel<SortmernaResult> = SORTMERNA(
            ch_reads,
            ch_sortmerna_fastas,
            ch_sortmerna_idx,
        )

        ch_results = ch_sortmerna
            .map { r ->
                record(
                    id:               r.id,
                    meta:             r.meta,
                    reads:            r.reads.isEmpty() ? null : r.reads,
                    sortmerna_log:    r.log,
                    ribodetector_log: null,
                    seqkit_stats:     null,
                    bowtie2_log:      null
                )
            }
    }
    else if (ribo_removal_tool == 'ribodetector') {
        // Run seqkit stats to determine average read length
        def ch_seqkit_stats: Channel<SeqkitStatsResult> = SEQKIT_STATS(ch_reads)

        // Join stats with reads and calculate read length for RiboDetector
        ch_reads_with_stats = ch_reads.join(ch_seqkit_stats, by: 'id')
        def ch_ribodetector: Channel<RibodetectorResult> = RIBODETECTOR(
            ch_reads_with_stats.map { r -> r + record(length: getReadLengthFromSeqkitStats(r.stats), gpu: tool_args.use_gpu_ribodetector ?: false) }
        )

        ch_results = ch_ribodetector
            .join(ch_seqkit_stats.map { r -> record(id: r.id, seqkit_stats: r.stats) }, by: 'id')
            .map { r ->
                record(
                    id:               r.id,
                    meta:             r.meta,
                    reads:            r.reads,
                    sortmerna_log:    null,
                    ribodetector_log: r.log,
                    seqkit_stats:     r.seqkit_stats,
                    bowtie2_log:      null
                )
            }
    }
    else {
        ch_bowtie2_idx = ch_bowtie2_index
        if (make_bowtie2_index) {
            // Process each rRNA file to add unique prefixes and convert U to T
            // This prevents duplicate sequence IDs in SAM header when combining databases
            ch_rrna_with_meta = ch_rrna_fastas.map { fasta_file ->
                record(id: fasta_file.baseName, meta: [id: fasta_file.baseName], fastx: fasta_file)
            }

            // Step 1: Add filename prefixes to sequence headers
            def ch_seqkit_prefixed: Channel<SeqkitReplaceResult> = SEQKIT_REPLACE(ch_rrna_with_meta, '')

            // Step 2: Convert U to T in sequences (RNA to DNA)
            ch_prefixed_fastas = ch_seqkit_prefixed.map { r ->
                record(id: "${r.meta.id}_dna", meta: [id: "${r.meta.id}_dna"], fastx: r.fastx)
            }
            def ch_seqkit_converted: Channel<SeqkitReplaceResult> = SEQKIT_REPLACE_U2T(ch_prefixed_fastas, '')

            // Collect processed files (already prefixed and U->T converted)
            def ch_combined_fasta = ch_seqkit_converted
                .map { r -> r.fastx }
                .collectFile(name: 'rrna_combined_dna.fasta', sort: { a, b -> a.name <=> b.name }, newLine: true)
                .collect()
                .map { fastas -> record(id: 'rrna_refs', meta: [id: 'rrna_refs'], fasta: fastas.toList().first() as Path) }

            def ch_bowtie2_built: Value<Bowtie2BuildResult> = BOWTIE2_BUILD(ch_combined_fasta)
            ch_bowtie2_idx = ch_bowtie2_built.map { built -> built.index }
            val_seqkit_prefixed = ch_seqkit_prefixed
                .collect()
                .map { built -> built.collect { r -> r.fastx }.toSorted { f -> f.name } }
            val_seqkit_converted = ch_seqkit_converted
                .collect()
                .map { built -> built.collect { r -> r.fastx }.toSorted { f -> f.name } }
            val_refs = ch_bowtie2_built
                .combine(val_seqkit_prefixed)
                .combine(val_seqkit_converted)
                .map { built, prefixed, converted ->
                    record(sortmerna_index: null, bowtie2_index: built.index, seqkit_prefixed: prefixed, seqkit_converted: converted)
                }
        }

        // For single-end reads: bowtie2's --un-gz works correctly
        // save_unaligned=true outputs unmapped reads directly
        def ch_bowtie2_se: Channel<Bowtie2AlignResult> = BOWTIE2_ALIGN(
            ch_reads.filter { r -> r.meta.single_end },
            ch_bowtie2_idx,
            null,             // No reference fasta needed
            true,             // save_unaligned - for single-end this works correctly
            false,            // sort_bam - not needed
        )

        // For paired-end reads: bowtie2's --un-conc-gz outputs pairs that didn't
        // align concordantly, which INCLUDES pairs where one mate aligned.
        // We need to filter via samtools to get pairs where BOTH mates are unmapped.
        def ch_bowtie2_pe: Channel<Bowtie2AlignResult> = BOWTIE2_ALIGN_PE(
            ch_reads.filter { r -> !r.meta.single_end },
            ch_bowtie2_idx,
            null,             // No reference fasta needed for BAM output
            false,            // save_unaligned - we'll extract from BAM instead
            false,            // sort_bam - not needed
        )

        // Filter BAM for read pairs where BOTH mates are unmapped (flag 12 = 4 + 8)
        // This removes any pair where at least one mate aligned to rRNA
        def ch_view: Channel<SamtoolsViewResult> = SAMTOOLS_VIEW_BOWTIE2(
            ch_bowtie2_pe.filter { r -> !r.raw_bams.isEmpty() }.map { r -> record(id: r.id, meta: r.meta, bam: r.raw_bams[0], bai: null) },
            null, // No reference fasta
            null, // No reference index
            null, // No qname file
            null, // No bed file
            ''    // No index format
        )
        // Note: samtools/view versions collected via topic
        ch_view_bam = ch_view.filter { r -> r.bam != null }

        // Convert filtered BAM back to paired FASTQ
        def ch_fastq_pe: Channel<SamtoolsFastqResult> = SAMTOOLS_FASTQ_BOWTIE2(
            ch_view_bam,
            false, // not interleaved
        )

        // Combine single-end and paired-end results
        ch_filtered_reads = ch_bowtie2_se
            .filter { r -> !r.unmapped.isEmpty() }
            .map { r -> record(id: r.id, reads: r.unmapped) }
            .mix(ch_fastq_pe.filter { r -> !r.reads.isEmpty() })

        ch_bowtie2_logs = ch_bowtie2_se
            .mix(ch_bowtie2_pe)
            .map { r ->
                record(
                    id:               r.id,
                    meta:             r.meta,
                    reads:            null,
                    sortmerna_log:    null,
                    ribodetector_log: null,
                    seqkit_stats:     null,
                    bowtie2_log:      r.bowtie2.log
                )
            }

        ch_results = ch_bowtie2_logs
            .join(ch_filtered_reads, by: 'id', remainder: true)
    }

    emit:
    samples: Channel<FastqRemoveRrna> = ch_results
    references: Value<RrnaReferences> = val_refs
}
