include { BBMAP_BBSPLIT                         } from '../../../modules/nf-core/bbmap/bbsplit'
include { FASTQC as FASTQC_FILTERED             } from '../../../modules/nf-core/fastqc'
include { CAT_FASTQ                             } from '../../../modules/nf-core/cat/fastq/main'
include { FQ_LINT                               } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_TRIMMING     } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_BBSPLIT      } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_RIBO_REMOVAL } from '../../../modules/nf-core/fq/lint/main'
include { FASTQ_REMOVE_RRNA                     } from '../fastq_remove_rrna'
include { FASTQ_SUBSAMPLE_FQ_SALMON             } from '../fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from '../fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP           } from '../fastq_fastqc_umitools_fastp'
include { FastqQcTrimFilterSetstrandedness      } from './types'

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

    // Calculate unstranded fragments (IU and U)
    // NOTE: this is here for completeness, but actually all fragments have a
    // strandedness (even if the overall library does not), so all these values
    // will be '0'. See
    // https://groups.google.com/g/sailfish-users/c/yxzBDv6NB6I
    def unstrandedKeys = ['IU', 'U', 'MU']
    def unstrandedFragments = unstrandedKeys.collect { key -> libCounts[key] ?: 0 }.sum()

    // The per-format keys (SF, ISF, MSF, OSF, ...) are not fragment-exclusive:
    // a fragment's candidate placements can increment both the sense and
    // antisense variant. `strand_mapping_bias` is salmon's fragment-exclusive
    // measure of the sense share of the resolved orientation instead.
    def numAssignedFragments = (libCounts['num_assigned_fragments'] ?: 0) as double
    def strandMappingBias = (libCounts['strand_mapping_bias'] ?: 0.0) as double
    def forwardFragments = strandMappingBias * numAssignedFragments
    def reverseFragments = (1 - strandMappingBias) * numAssignedFragments

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
    // Input channels
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: /path/to/genome.fasta
    ch_transcript_fasta  // channel: /path/to/transcript.fasta
    ch_gtf               // channel: /path/to/genome.gtf
    ch_salmon_index      // channel: /path/to/salmon/index/ (optional)
    ch_sortmerna_index   // channel: /path/to/sortmerna/index/ (optional)
    ch_bowtie2_index     // channel: /path/to/bowtie2/index/ (optional)
    ch_bbsplit_index     // channel: /path/to/bbsplit/index/ (optional)
    ch_rrna_fastas       // channel: one or more fasta files containing rrna sequences to be passed to SortMeRNA/Bowtie2 (optional)

    // Skip options
    skip_bbsplit         // boolean: Skip BBSplit for removal of non-reference genome reads.
    skip_fastqc          // boolean: true/false
    skip_trimming        // boolean: true/false
    skip_umi_extract     // boolean: true/false
    skip_linting         // boolean: true/false

    // Index generation
    make_salmon_index    // boolean: Whether to create salmon index before running salmon quant
    make_sortmerna_index // boolean: Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index   // boolean: Whether to create a bowtie2 index before running bowtie2

    // Trimming options
    trimmer              // string (enum): 'fastp' or 'trimgalore'
    min_trimmed_reads    // integer: > 0
    save_trimmed         // boolean: true/false
    fastp_merge          // boolean: true/false: whether to stitch paired end reads together in FASTP output

    // rRNA removal options
    remove_ribo_rna           // boolean: true/false: whether to remove rRNA
    ribo_removal_tool         // string (enum): 'sortmerna', 'ribodetector', or 'bowtie2'

    // UMI options
    with_umi             // boolean: true/false: Enable UMI-based read deduplication.
    umi_discard_read     // integer: 0, 1 or 2

    // Merging options
    save_merged_fastq    // boolean: true/false: Save merged FastQ files even for single-library samples

    // Strandedness thresholds
    stranded_threshold   // float: The fraction of stranded reads that must be assigned to a strandedness for confident assignment. Must be at least 0.5
    unstranded_threshold // float: The difference in fraction of stranded reads assigned to 'forward' and 'reverse' below which a sample is classified as 'unstranded'

    main:

    ch_filtered_reads = channel.empty()
    ch_trim_read_count = channel.empty()
    ch_trim_reads_merged = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_lint_log_raw = channel.empty()
    ch_lint_log_trimmed = channel.empty()
    ch_lint_log_bbsplit = channel.empty()
    ch_lint_log_ribo = channel.empty()

    // Individual output channels for workflow outputs
    ch_fastqc_raw_html    = channel.empty()
    ch_fastqc_raw_zip     = channel.empty()
    ch_fastqc_trim_html   = channel.empty()
    ch_fastqc_trim_zip    = channel.empty()
    ch_trim_html          = channel.empty()
    ch_trim_zip           = channel.empty()
    ch_trim_log           = channel.empty()
    ch_trim_json          = channel.empty()
    ch_trim_unpaired      = channel.empty()
    ch_umi_log            = channel.empty()
    ch_umi_reads          = channel.empty()
    ch_bbsplit_stats      = channel.empty()
    ch_sortmerna_log      = channel.empty()
    ch_ribodetector_log   = channel.empty()
    ch_seqkit_stats       = channel.empty()
    ch_bowtie2_log        = channel.empty()
    ch_sortmerna_index_built = channel.empty()
    ch_salmon_index_built    = channel.empty()
    ch_seqkit_prefixed    = channel.empty()
    ch_seqkit_converted   = channel.empty()
    ch_fastqc_filtered_html = channel.empty()
    ch_fastqc_filtered_zip  = channel.empty()

    ch_reads
        .branch { meta, fastqs ->
            single: fastqs.size() == 1 && fastqs.flatten()[0].name.endsWith('.gz') && !save_merged_fastq
            return [meta, fastqs.flatten()]
            multiple: true
            return [meta, fastqs.flatten()]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(
        ch_fastq.multiple
    ).reads.mix(ch_fastq.single).set { ch_filtered_reads }

    //
    // MODULE: Lint FastQ files
    //

    if (!skip_linting) {
        FQ_LINT(
            ch_filtered_reads
        )
        ch_lint_log_raw = FQ_LINT.out.lint
        ch_filtered_reads = ch_filtered_reads.join(FQ_LINT.out.lint.map { meta, _lint -> meta })
    }

    ch_reads_cat = ch_filtered_reads

    // Each stage that runs joins its outputs onto this per-sample skeleton;
    // fields for skipped stages are left unset and become null in the final
    // record, so no join is ever made against an empty channel. Stages after
    // the min_trimmed_reads filter join with remainder so that samples
    // failing it keep their record.
    ch_results = ch_reads_cat.map { meta, reads -> [meta.id, [meta: meta, reads_cat: [reads].flatten()]] }

    if (!skip_linting) {
        ch_results = ch_results
            .join(ch_lint_log_raw.map { meta, lint -> [meta.id, lint] }, by: [0])
            .map { id, fields, lint -> [id, fields + [lint_raw: lint]] }
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

        // Capture individual outputs for workflow outputs
        ch_fastqc_raw_html = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_html
        ch_fastqc_raw_zip  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_html = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_html
        ch_fastqc_trim_zip  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_json       = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_json
        ch_trim_log        = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_unpaired   = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_unpaired
        ch_umi_log         = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.umi_log
        ch_umi_reads       = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.umi_reads

        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_json)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.umi_log)
            .mix(ch_multiqc_files)

        ch_results = ch_results
            .join(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.results.map { r -> [r.id, r] }, by: [0])
            .map { id, fields, r ->
                def trim = r.trim
                    ? record(html: null, log: r.trim.log, json: r.trim.json, unpaired: r.trim.unpaired, reads_fail: null, reads_merged: null)
                    : null
                [id, fields + [
                    fastqc_raw_html:   r.fastqc?.raw_html,
                    fastqc_raw_zip:    r.fastqc?.raw_zip,
                    fastqc_trim_html:  r.trim?.html,
                    fastqc_trim_zip:   r.trim?.zip,
                    trim:              trim,
                    umi:               r.umi ? record(log: r.umi.log, reads: r.umi.reads) : null,
                    num_trimmed_reads: r.num_trimmed_reads != null ? r.num_trimmed_reads as Long : null
                ]]
            }
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP(
            ch_filtered_reads.map { meta, reads -> tuple(meta, reads, []) }, // Add empty adapter sequence
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            save_trimmed,
            fastp_merge,
            min_trimmed_reads,
        )
        ch_filtered_reads = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_trim_reads_merged = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_reads_merged

        // Capture individual outputs for workflow outputs
        ch_fastqc_raw_html  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_html
        ch_fastqc_raw_zip   = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_html = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_html
        ch_fastqc_trim_zip  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_json        = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_html        = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_html
        ch_trim_log         = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_log
        ch_umi_log          = FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_log
        ch_umi_reads        = FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_reads

        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_log)
            .mix(ch_multiqc_files)

        ch_results = ch_results
            .join(FASTQ_FASTQC_UMITOOLS_FASTP.out.results.map { r -> [r.id, r] }, by: [0])
            .map { id, fields, r ->
                def trim = r.trim
                    ? record(html: [r.trim.html], log: [r.trim.log], json: [r.trim.json], unpaired: null, reads_fail: r.trim.reads_fail, reads_merged: r.trim.reads_merged)
                    : null
                [id, fields + [
                    fastqc_raw_html:   r.fastqc?.raw_html,
                    fastqc_raw_zip:    r.fastqc?.raw_zip,
                    fastqc_trim_html:  r.fastqc?.trim_html,
                    fastqc_trim_zip:   r.fastqc?.trim_zip,
                    trim:              trim,
                    umi:               r.umi ? record(log: r.umi.log, reads: r.umi.reads) : null,
                    num_trimmed_reads: r.num_trimmed_reads
                ]]
            }
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
        ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').map { file -> [[:], file] }
    )

    if ((!skip_linting) && (!skip_trimming)) {
        FQ_LINT_AFTER_TRIMMING(
            ch_filtered_reads
        )
        ch_lint_log_trimmed = FQ_LINT_AFTER_TRIMMING.out.lint
        ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_TRIMMING.out.lint.map { meta, _lint -> meta })

        ch_results = ch_results
            .join(ch_lint_log_trimmed.map { meta, lint -> [meta.id, lint] }, by: [0], remainder: true)
            .map { id, fields, lint -> [id, fields + [lint_trimmed: lint]] }
    }

    ch_reads_trimmed = ch_filtered_reads

    ch_results = ch_results
        .join(ch_reads_trimmed.map { meta, reads -> [meta.id, [reads].flatten()] }, by: [0], remainder: true)
        .map { id, fields, reads -> [id, fields + [reads_trimmed: reads]] }

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

        ch_bbsplit_stats = BBMAP_BBSPLIT.out.stats
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBSPLIT.out.stats)

        // all_fastq includes the primary reads too; keep only the rest for publishing.
        ch_bbsplit_other_reads = BBMAP_BBSPLIT.out.all_fastq
            .join(BBMAP_BBSPLIT.out.primary_fastq)
            .map { meta, all, primary ->
                def primary_names = [primary].flatten()*.name
                [meta, [all].flatten().findAll { f -> !(f.name in primary_names) }]
            }

        ch_results = ch_results
            .join(ch_bbsplit_stats.map { meta, stats -> [meta.id, stats] }, by: [0], remainder: true)
            .join(ch_bbsplit_other_reads.map { meta, other -> [meta.id, other] }, by: [0], remainder: true)
            .join(BBMAP_BBSPLIT.out.primary_fastq.map { meta, primary -> [meta.id, [primary].flatten()] }, by: [0], remainder: true)
            .map { id, fields, stats, other_reads, primary_reads -> [id, fields + [bbsplit: stats != null ? record(stats: stats, primary_reads: primary_reads, other_genome_reads: other_reads) : null]] }

        if (!skip_linting) {
            FQ_LINT_AFTER_BBSPLIT(
                ch_filtered_reads
            )
            ch_lint_log_bbsplit = FQ_LINT_AFTER_BBSPLIT.out.lint
            ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_BBSPLIT.out.lint.map { meta, _lint -> meta })

            ch_results = ch_results
                .join(ch_lint_log_bbsplit.map { meta, lint -> [meta.id, lint] }, by: [0], remainder: true)
                .map { id, fields, lint -> [id, fields + [lint_bbsplit: lint]] }
        }
    }

    //
    // SUBWORKFLOW: Remove ribosomal RNA reads
    //
    if (remove_ribo_rna) {
        FASTQ_REMOVE_RRNA(
            ch_filtered_reads,
            ch_rrna_fastas,
            ch_sortmerna_index,
            ch_bowtie2_index,
            ribo_removal_tool,
            make_sortmerna_index,
            make_bowtie2_index,
        )

        ch_filtered_reads = FASTQ_REMOVE_RRNA.out.reads
        ch_sortmerna_log    = FASTQ_REMOVE_RRNA.out.sortmerna_log
        ch_ribodetector_log = FASTQ_REMOVE_RRNA.out.ribodetector_log
        ch_seqkit_stats     = FASTQ_REMOVE_RRNA.out.seqkit_stats
        ch_bowtie2_log      = FASTQ_REMOVE_RRNA.out.bowtie2_log
        ch_bowtie2_index    = FASTQ_REMOVE_RRNA.out.bowtie2_index
        ch_sortmerna_index_built = FASTQ_REMOVE_RRNA.out.sortmerna_index
        ch_seqkit_prefixed  = FASTQ_REMOVE_RRNA.out.seqkit_prefixed
        ch_seqkit_converted = FASTQ_REMOVE_RRNA.out.seqkit_converted
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_REMOVE_RRNA.out.multiqc_files)

        ch_results = ch_results
            .join(FASTQ_REMOVE_RRNA.out.results.map { r -> [r.id, r] }, by: [0], remainder: true)
            .map { id, fields, r -> [id, fields + [rrna: r]] }

        if (!skip_linting) {
            FQ_LINT_AFTER_RIBO_REMOVAL(
                ch_filtered_reads
            )
            ch_lint_log_ribo = FQ_LINT_AFTER_RIBO_REMOVAL.out.lint
            ch_filtered_reads = ch_filtered_reads.join(FQ_LINT_AFTER_RIBO_REMOVAL.out.lint.map { meta, _lint -> meta })

            ch_results = ch_results
                .join(ch_lint_log_ribo.map { meta, lint -> [meta.id, lint] }, by: [0], remainder: true)
                .map { id, fields, lint -> [id, fields + [lint_ribo: lint]] }
        }
    }

    //
    // MODULE: Run FastQC on filtered reads (after BBSplit and/or rRNA removal)
    //
    if (!skip_fastqc && (!skip_bbsplit || remove_ribo_rna)) {
        FASTQC_FILTERED(
            ch_filtered_reads
        )
        ch_fastqc_filtered_html = FASTQC_FILTERED.out.html
        ch_fastqc_filtered_zip  = FASTQC_FILTERED.out.zip
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_FILTERED.out.zip)

        ch_results = ch_results
            .join(FASTQC_FILTERED.out.html.map { meta, html -> [meta.id, [html].flatten()] }, by: [0], remainder: true)
            .join(FASTQC_FILTERED.out.zip.map { meta, zip -> [meta.id, [zip].flatten()] }, by: [0], remainder: true)
            .map { id, fields, html, zip -> [id, fields + [fastqc_filtered_html: html, fastqc_filtered_zip: zip]] }
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
        .map { fasta -> [fasta] }
        .combine(ch_strand_fastq.auto_strand)
        .map { items -> items[0] }
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

    ch_salmon_index_built = FASTQ_SUBSAMPLE_FQ_SALMON.out.index_built

    FASTQ_SUBSAMPLE_FQ_SALMON.out.lib_format_counts
        .join(ch_strand_fastq.auto_strand, remainder: true)
        .map { meta, json, reads ->
            if (json == null) {
                error("Salmon failed to produce lib_format_counts for sample '${meta.id}' " +
                    "which was set to 'auto' strandedness. Check that the Salmon " +
                    "index matches your input reads, or set strandedness explicitly in the samplesheet.")
            }
            def salmon_strand_analysis = getSalmonInferredStrandedness(json, stranded_threshold, unstranded_threshold)
            def strandedness = salmon_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                strandedness = 'unstranded'
            }
            return [meta + [strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis], reads]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    ch_results = ch_results
        .join(ch_strand_inferred_fastq.map { meta, reads -> [meta.id, [meta: meta, reads: [reads].flatten()]] }, by: [0], remainder: true)
        .map { id, fields, final_reads ->
            def has_fastqc = [
                fields.fastqc_raw_html, fields.fastqc_raw_zip,
                fields.fastqc_trim_html, fields.fastqc_trim_zip,
                fields.fastqc_filtered_html, fields.fastqc_filtered_zip
            ].any { f -> f != null }
            def fastqc = has_fastqc
                ? record(
                    raw_html:      fields.fastqc_raw_html,
                    raw_zip:       fields.fastqc_raw_zip,
                    trim_html:     fields.fastqc_trim_html,
                    trim_zip:      fields.fastqc_trim_zip,
                    filtered_html: fields.fastqc_filtered_html,
                    filtered_zip:  fields.fastqc_filtered_zip
                )
                : null
            def lint = skip_linting
                ? null
                : record(
                    raw:     fields.lint_raw,
                    trimmed: fields.lint_trimmed,
                    bbsplit: fields.lint_bbsplit,
                    ribo:    fields.lint_ribo
                )
            record(
                id:                id,
                meta:              final_reads?.meta ?: fields.meta,
                reads:             final_reads?.reads,
                reads_cat:         fields.reads_cat,
                reads_trimmed:     fields.reads_trimmed,
                num_trimmed_reads: fields.num_trimmed_reads,
                fastqc:            fastqc,
                trim:              fields.trim,
                umi:               fields.umi,
                bbsplit:           fields.bbsplit,
                lint:              lint,
                rrna:              fields.rrna
            )
        }

    // `remainder: true` needed because every contributor is gated behind
    // a `skip_*` / `filter_*` option.
    ch_per_sample_mqc_bundle = ch_fastqc_raw_zip
        .join(ch_fastqc_trim_zip,      remainder: true)
        .join(ch_trim_log,             remainder: true)
        .join(ch_trim_json,            remainder: true)
        .join(ch_umi_log,              remainder: true)
        .join(ch_bbsplit_stats,        remainder: true)
        .join(ch_sortmerna_log,        remainder: true)
        .join(ch_ribodetector_log,     remainder: true)
        .join(ch_seqkit_stats,         remainder: true)
        .join(ch_bowtie2_log,          remainder: true)
        .join(ch_fastqc_filtered_zip,  remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    emit:
    reads             = ch_strand_inferred_fastq
    reads_cat         = ch_reads_cat
    reads_trimmed     = ch_reads_trimmed
    trim_read_count   = ch_trim_read_count
    trim_reads_merged = ch_trim_reads_merged
    multiqc_files     = ch_multiqc_files.transpose()

    // Individual outputs for workflow outputs
    lint_log_raw     = ch_lint_log_raw
    lint_log_trimmed = ch_lint_log_trimmed
    lint_log_bbsplit = ch_lint_log_bbsplit
    lint_log_ribo    = ch_lint_log_ribo
    fastqc_raw_html  = ch_fastqc_raw_html
    fastqc_raw_zip   = ch_fastqc_raw_zip
    fastqc_trim_html = ch_fastqc_trim_html
    fastqc_trim_zip  = ch_fastqc_trim_zip
    trim_html        = ch_trim_html
    trim_zip         = ch_trim_zip
    trim_log         = ch_trim_log
    trim_json        = ch_trim_json
    trim_unpaired    = ch_trim_unpaired
    umi_log          = ch_umi_log
    umi_reads        = ch_umi_reads
    bbsplit_stats    = ch_bbsplit_stats
    sortmerna_log    = ch_sortmerna_log
    ribodetector_log = ch_ribodetector_log
    seqkit_stats     = ch_seqkit_stats
    bowtie2_log      = ch_bowtie2_log
    bowtie2_index    = ch_bowtie2_index
    sortmerna_index_built = ch_sortmerna_index_built
    salmon_index_built    = ch_salmon_index_built
    fastqc_filtered_html = ch_fastqc_filtered_html
    fastqc_filtered_zip  = ch_fastqc_filtered_zip
    seqkit_prefixed  = ch_seqkit_prefixed
    seqkit_converted = ch_seqkit_converted
    per_sample_mqc_bundle = ch_per_sample_mqc_bundle // channel: [ val(meta), list(files) ]
    results          = ch_results // channel: FastqQcTrimFilterSetstrandedness
}
