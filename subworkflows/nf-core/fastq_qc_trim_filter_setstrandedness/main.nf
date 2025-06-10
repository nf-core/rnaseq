include { BBMAP_BBSPLIT                      } from '../../../modules/nf-core/bbmap/bbsplit'
include { CAT_FASTQ                          } from '../../../modules/nf-core/cat/fastq/main'
include { SORTMERNA                          } from '../../../modules/nf-core/sortmerna/main'
include { SORTMERNA as SORTMERNA_INDEX       } from '../../../modules/nf-core/sortmerna/main'
include { FQ_LINT                            } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_TRIMMING  } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_BBSPLIT   } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_SORTMERNA } from '../../../modules/nf-core/fq/lint/main'

include { FASTQ_SUBSAMPLE_FQ_SALMON          } from '../fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE   } from '../fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP        } from '../fastq_fastqc_umitools_fastp'

//
// Function to determine library type by comparing type counts.
//

//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments

    def strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)

        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        }
        else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        }
        else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }

    return [
        inferred_strandedness: strandedness,
        forwardFragments: (forwardFragments / (totalFragments as double)) * 100,
        reverseFragments: (reverseFragments / (totalFragments as double)) * 100,
        unstrandedFragments: (unstrandedFragments / (totalFragments as double)) * 100,
    ]
}

//
// Function that parses Salmon quant 'lib_format_counts.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    // Parse the JSON content of the file
    def libCounts = new groovy.json.JsonSlurper().parseText(json_file.text)

    // Calculate the counts for forward and reverse strand fragments
    def forwardKeys = ['SF', 'ISF', 'MSF', 'OSF']
    def reverseKeys = ['SR', 'ISR', 'MSR', 'OSR']

    // Calculate unstranded fragments (IU and U)
    // NOTE: this is here for completeness, but actually all fragments have a
    // strandedness (even if the overall library does not), so all these values
    // will be '0'. See
    // https://groups.google.com/g/sailfish-users/c/yxzBDv6NB6I
    def unstrandedKeys = ['IU', 'U', 'MU']

    def forwardFragments = forwardKeys.collect { libCounts[it] ?: 0 }.sum()
    def reverseFragments = reverseKeys.collect { libCounts[it] ?: 0 }.sum()
    def unstrandedFragments = unstrandedKeys.collect { libCounts[it] ?: 0 }.sum()

    // Use shared calculation function to determine strandedness
    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}

workflow FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS {
    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: /path/to/genome.fasta
    ch_transcript_fasta  // channel: /path/to/transcript.fasta
    ch_gtf               // channel: /path/to/genome.gtf
    ch_salmon_index      // channel: /path/to/salmon/index/ (optional)
    ch_sortmerna_index   // channel: /path/to/sortmerna/index/ (optional)
    ch_bbsplit_index     // channel: /path/to/bbsplit/index/ (optional)
    ch_rrna_fastas       // channel: one or more fasta files containing rrna sequences to be passed to SortMeRNA (optional)
    skip_bbsplit         // boolean: Skip BBSplit for removal of non-reference genome reads.
    skip_fastqc          // boolean: true/false
    skip_trimming        // boolean: true/false
    skip_umi_extract     // boolean: true/false
    make_salmon_index    // boolean: Whether to create salmon index before running salmon quant
    make_sortmerna_index // boolean: Whether to create a sortmerna index before running sortmerna
    trimmer              // string (enum): 'fastp' or 'trimgalore'
    min_trimmed_reads    // integer: > 0
    save_trimmed         // boolean: true/false
    remove_ribo_rna      // boolean: true/false: whether to run sortmerna to remove rrnas
    with_umi             // boolean: true/false: Enable UMI-based read deduplication.
    umi_discard_read     // integer: 0, 1 or 2
    stranded_threshold   // float: The fraction of stranded reads that must be assigned to a strandedness for confident assignment. Must be at least 0.5
    unstranded_threshold // float: The difference in fraction of stranded reads assigned to 'forward' and 'reverse' below which a sample is classified as 'unstranded'
    skip_linting         // boolean: true/false
    fastp_merge          // boolean: true/false: whether to stitch paired end reads together in FASTP output

    main:

    ch_versions = Channel.empty()
    ch_filtered_reads = Channel.empty()
    ch_trim_read_count = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_lint_log = Channel.empty()

    ch_reads
        .branch { meta, fastqs ->
            single: fastqs.size() == 1
            return [meta, fastqs.flatten()]
            multiple: fastqs.size() > 1
            return [meta, fastqs.flatten()]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(
        ch_fastq.multiple
    ).reads.mix(ch_fastq.single).set { ch_filtered_reads }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Lint FastQ files
    //

    if (!skip_linting) {
        FQ_LINT(
            ch_filtered_reads
        )
        ch_versions = ch_versions.mix(FQ_LINT.out.versions.first())
        ch_lint_log = ch_lint_log.mix(FQ_LINT.out.lint)
        ch_filtered_reads = ch_filtered_reads.join(FQ_LINT.out.lint.map { it[0] })
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    if (trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
            ch_filtered_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            skip_trimming,
            umi_discard_read,
            min_trimmed_reads,
        )
        ch_filtered_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .mix(ch_multiqc_files)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP(
            ch_filtered_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            [],
            save_trimmed,
            fastp_merge,
            min_trimmed_reads,
        )
        ch_filtered_reads = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
            .mix(ch_multiqc_files)
    }

    def pass_trimmed_reads = [:]

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //

    ch_trim_read_count
        .map { meta, num_reads ->
            pass_trimmed_reads[meta.id] = true
            if (num_reads <= min_trimmed_reads.toFloat()) {
                pass_trimmed_reads[meta.id] = false
                return ["${meta.id}\t${num_reads}"]
            }
        }
        .collect()
        .map { tsv_data ->
            def header = ["Sample", "Reads after trimming"]
            multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').map { [[:], it] }
    )

    if ((!skip_linting) && (!skip_trimming)) {
        FQ_LINT_AFTER_TRIMMING(
            ch_filtered_reads
        )
        ch_lint_log = ch_lint_log.mix(FQ_LINT_AFTER_TRIMMING.out.lint)
        ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_TRIMMING.out.lint.map { it[0] })
    }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!skip_bbsplit) {
        BBMAP_BBSPLIT(
            ch_filtered_reads,
            ch_bbsplit_index,
            [],
            [[], []],
            false,
        )

        BBMAP_BBSPLIT.out.primary_fastq.set { ch_filtered_reads }

        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())

        if (!skip_linting) {
            FQ_LINT_AFTER_BBSPLIT(
                ch_filtered_reads
            )
            ch_lint_log = ch_lint_log.mix(FQ_LINT_AFTER_BBSPLIT.out.lint)
            ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_BBSPLIT.out.lint.map { it[0] })
        }
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    if (remove_ribo_rna) {
        ch_sortmerna_fastas = ch_rrna_fastas
            .collect()
            .map { ['rrna_refs', it] }

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

        SORTMERNA.out.reads.set { ch_filtered_reads }

        ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log)

        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())

        if (!skip_linting) {
            FQ_LINT_AFTER_SORTMERNA(
                ch_filtered_reads
            )
            ch_lint_log = ch_lint_log.mix(FQ_LINT_AFTER_SORTMERNA.out.lint)
            ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_SORTMERNA.out.lint.map { it[0] })
        }
    }

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_filtered_reads
        .branch { meta, fastq ->
            auto_strand: meta.strandedness == 'auto'
            return [meta, fastq]
            known_strand: meta.strandedness != 'auto'
            return [meta, fastq]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created

    ch_fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON(
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index,
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON.out.lib_format_counts
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            def salmon_strand_analysis = getSalmonInferredStrandedness(json, stranded_threshold, unstranded_threshold)
            def strandedness = salmon_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                strandedness = 'unstranded'
            }
            return [meta + [strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis], reads]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    emit:
    lint_log        = ch_lint_log
    reads           = ch_strand_inferred_fastq
    trim_read_count = ch_trim_read_count
    multiqc_files   = ch_multiqc_files.transpose().map { it[1] }
    versions        = ch_versions // channel: [ versions.yml ]
}
