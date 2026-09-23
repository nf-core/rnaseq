//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { samplesheetToList      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'
include { workflowVersionToYAML  } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../utils_nfcore_rnaseq_pipeline'


workflow MULTIQC_RNASEQ {

    take:
    ch_multiqc_files           // channel: [ val(meta), path(file) ]       - flat, contributor outputs
    ch_per_sample_bundle_raw   // channel: [ id, meta, f1, f2, ... ]       - per-sample, grown by `.join(..., remainder: true)` at each subworkflow aggregation site
    ch_strand_data             // channel: [ val(meta), provided, status, salmon, rseqc ] - per-sample strand classification, used for the Strandedness checks section
    ch_trim_read_count         // channel: [ val(meta), val(num_reads) ]   - for fail_trimmed section
    ch_percent_mapped_pass     // channel: [ id, percent_mapped, pass ]    - for fail_mapped section
    aligner_display_name       // string: display name of the aligner used for the percent_mapped metric, e.g. 'STAR uniquely mapped reads' or 'Bowtie2 overall alignment rate'
    ch_fastq                   // channel: [ val(meta), [ reads ] ]
    ch_collated_versions       // channel: path(versions yaml)
    samplesheet_path           // path: pipeline input samplesheet
    samplesheet_schema         // path: samplesheet JSON schema
    mqc_default_config         // path: pipeline-bundled MultiQC config
    mqc_custom_config          // path (or []): optional user MultiQC config
    mqc_logo                   // path (or []): optional custom logo
    methods_description_yml    // path: methods-description YAML template
    strand_summary_asset       // path: strand_check_summary YAML custom-content template
    strand_composition_asset   // path: strand_check_composition YAML custom-content template
    sample_status_header       // path: MultiQC custom content header for fail_* tables
    min_trimmed_reads          // integer: threshold for fail_trimmed classification
    skip_quantification_merge  // boolean

    main:

    //
    // fail_* custom-content TSVs. Each sample either contributes a
    // single fail row or an empty placeholder so the downstream
    // per-sample `.join(..., remainder: true)` can close
    // progressively. The anchor is derived from the bundle itself so
    // every bundle sample has a match on every fail_* stream.
    //
    // `status_header_lines` tracks the header row count so editing
    // `sample_status_header.txt` doesn't silently mis-skip the merged
    // aggregate's concatenation.
    //
    def status_header_lines = sample_status_header.readLines().size() + 1  // parent header + one column row
    ch_sample_anchor_by_id  = ch_per_sample_bundle_raw.map { row -> [row[0], row[1]] }

    ch_fail_trimmed_fail_by_id = ch_trim_read_count
        .filter { _meta, n -> n <= min_trimmed_reads.toFloat() }
        .collectFile { meta, n ->
            [
                "${meta.id}_fail_trimmed_samples_mqc.tsv",
                "Sample\tReads after trimming\n${meta.id}\t${n}\n",
            ]
        }
        .map { f -> [f.baseName.replace('_fail_trimmed_samples_mqc', ''), f] }

    ch_fail_trimmed_all = ch_sample_anchor_by_id
        .join(ch_fail_trimmed_fail_by_id, remainder: true)
        .map { _id, meta, f -> [meta, f ?: []] }

    ch_fail_trimmed_merged = ch_fail_trimmed_all
        .map { _meta, f -> f }
        .flatten()
        .collectFile(name: 'fail_trimmed_samples_mqc.tsv', keepHeader: true)
        .map { f -> [[:], f] }

    ch_fail_mapped_fail_by_id = ch_percent_mapped_pass
        .filter { _id, _pm, pass -> pass != null && !pass }
        .collectFile { id, percent_mapped, _pass ->
            [
                "${id}_fail_mapped_samples_mqc.tsv",
                sample_status_header.text + "Sample\t${aligner_display_name} (%)\n${id}\t${percent_mapped}\n",
            ]
        }
        .map { f -> [f.baseName.replace('_fail_mapped_samples_mqc', ''), f] }

    ch_fail_mapped_all = ch_sample_anchor_by_id
        .join(ch_fail_mapped_fail_by_id, remainder: true)
        .map { _id, meta, f -> [meta, f ?: []] }

    ch_fail_mapped_merged = ch_fail_mapped_all
        .map { _meta, f -> f }
        .flatten()
        .collectFile(name: 'fail_mapped_samples_mqc.tsv', keepHeader: true, skip: status_header_lines)
        .map { f -> [[:], f] }

    //
    // Strandedness checks custom-content section. Two MultiQC
    // subsections (summary table + stacked composition bargraph) are
    // rendered from the same per-sample tuple, with header / pconfig
    // / colour config in the bundled YAML templates. The composition
    // section inherits `parent_*` from the summary section so the
    // description lives in one place.
    //
    def strand_summary_static     = loadMultiqcAsset(strand_summary_asset)
    def strand_composition_static = loadMultiqcAsset(strand_composition_asset) + strand_summary_static.subMap(['parent_id', 'parent_name', 'parent_description'])

    // Per-run table_sample_merge config: only PE samples from the
    // samplesheet get their _1 / _2 rows grouped in the General Stats
    // table.
    ch_mqc_dynamic_config = channel.of(multiqcSampleMergeYaml(samplesheet_path, samplesheet_schema))
        .collectFile(name: 'multiqc_sample_merge.yml')

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel
        .value(paramsSummaryMultiqc(paramsSummaryMap(workflow, parameters_schema: 'nextflow_schema.json')))
        .collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description = channel
        .value(methodsDescriptionText(methods_description_yml))
        .collectFile(name: 'methods_description_mqc.yaml')

    //
    // Two execution modes for MULTIQC:
    //   - merged (default): one report covers the whole run.
    //   - per-sample (--skip_quantification_merge): one report per
    //     sample; workflow-level versions are replaced with a
    //     pipeline-identity manifest so the report doesn't wait on
    //     the global versions topic.
    //
    // Each branch ends with a tuple matching the MULTIQC input
    // contract (id, files, configs, logo, replace_names, extra); the
    // closure below builds it so the branches stay focused on file
    // assembly.
    //
    def buildMultiqcInputTuple = { id, files, dynamic_config, replace_names = [] ->
        [
            [id: id],
            files,
            [mqc_default_config, dynamic_config, mqc_custom_config].findAll { cfg -> cfg },
            mqc_logo,
            replace_names,
            [],
        ]
    }

    if (skip_quantification_merge) {
        ch_strand_summary_by_id = ch_strand_data
            .collectFile { row ->
                [
                    "${row[0].id}_strand_check_summary_mqc.json",
                    strandCheckSummaryYaml(strand_summary_static, [row]),
                ]
            }
            .map { f -> [f.baseName.replace('_strand_check_summary_mqc', ''), f] }

        ch_strand_composition_by_id = ch_strand_data
            .collectFile { row ->
                [
                    "${row[0].id}_strand_check_composition_mqc.json",
                    strandCheckCompositionYaml(strand_composition_static, [row]),
                ]
            }
            .map { f -> [f.baseName.replace('_strand_check_composition_mqc', ''), f] }

        // Collapse the raw bundle with every per-sample contributor,
        // one `.join(remainder: true)` per stream. Each sample becomes
        // `[meta, [files]]`; missing streams show up as null entries
        // that are filtered out before MULTIQC sees them.
        ch_per_sample_bundle = ch_per_sample_bundle_raw
            .join(ch_fail_trimmed_all.map { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_fail_mapped_all.map  { meta, f -> [meta.id, f] }, remainder: true)
            .join(ch_strand_summary_by_id,     remainder: true)
            .join(ch_strand_composition_by_id, remainder: true)
            .map { row ->
                [
                    row[1],
                    row.drop(2)
                        .findAll { entry -> entry != null }
                        .collectMany { entry -> (entry instanceof List) ? entry : [entry] },
                ]
            }

        ch_manifest_versions = channel.value(workflowVersionToYAML())
            .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml')

        ch_static_globals = ch_workflow_summary
            .mix(ch_methods_description)
            .mix(ch_manifest_versions)
            .collect()

        ch_global_files = ch_fail_trimmed_merged
            .mix(ch_fail_mapped_merged)
            .map { _meta, f -> f }
            .collect()
            .ifEmpty([])

        ch_multiqc_input = ch_per_sample_bundle
            .combine(ch_static_globals.toList())
            .combine(ch_global_files.toList())
            .combine(ch_mqc_dynamic_config)
            .map { meta, sample_files, static_globals, run_globals, dyn ->
                // No replace_names: each per-sample report contains one sample.
                buildMultiqcInputTuple.call(
                    meta.id,
                    sample_files + (static_globals ?: []) + (run_globals ?: []),
                    dyn,
                )
            }
    } else {
        // `.collect(flat: false)` is silent on an empty channel, so
        // zero strand rows -> no *_mqc.json emission -> MultiQC drops
        // the section cleanly.
        ch_strand_rows = ch_strand_data.collect(flat: false)

        ch_strand_summary_merged = ch_strand_rows
            .map { rows -> strandCheckSummaryYaml(strand_summary_static, rows) }
            .collectFile(name: 'strand_check_summary_mqc.json')
            .map { f -> [[:], f] }

        ch_strand_composition_merged = ch_strand_rows
            .map { rows -> strandCheckCompositionYaml(strand_composition_static, rows) }
            .collectFile(name: 'strand_check_composition_mqc.json')
            .map { f -> [[:], f] }

        // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
        ch_name_replacements = multiqcNameReplacements(ch_fastq)

        // `multiqc_report` is a sentinel meta.id used by
        // conf/modules/multiqc.config to pick the merged output path.
        ch_multiqc_files_merged = ch_multiqc_files
            .mix(ch_fail_trimmed_merged)
            .mix(ch_fail_mapped_merged)
            .mix(ch_strand_summary_merged)
            .mix(ch_strand_composition_merged)
            .mix(ch_workflow_summary.mix(ch_collated_versions).mix(ch_methods_description).map { f -> [[:], f] })

        ch_multiqc_input = ch_multiqc_files_merged
            .map { _meta, f -> f }
            .collect()
            .map { files -> [files] }
            .combine(ch_name_replacements.ifEmpty([]).toList())
            .combine(ch_mqc_dynamic_config)
            .map { files, replace_names, dyn ->
                buildMultiqcInputTuple.call('multiqc_report', files, dyn, replace_names ?: [])
            }
    }

    MULTIQC(ch_multiqc_input)

    emit:
    report = MULTIQC.out.report.map { _meta, report -> report }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MultiQC `--replace-names` file: map each FASTQ simpleName to
// '<id>_1' / '<id>_2' (or '<id>' for SE), skipping cases where the
// simpleName already equals the sample ID (see #1341 / #1659).
//
def multiqcNameReplacements(ch_fastq) {
    return ch_fastq
        .map { meta, reads ->
            def paired   = reads[0][1] as boolean
            def suffixes = paired ? ['_1', '_2'] : ['']
            def mappings = []

            def fastq1_simplename = file(reads[0][0]).simpleName
            if (fastq1_simplename != meta.id) {
                mappings << [fastq1_simplename, "${meta.id}${suffixes[0]}"]
                if (paired) {
                    mappings << [file(reads[0][1]).simpleName, "${meta.id}${suffixes[1]}"]
                }
            }

            return mappings.collect { mapping -> mapping.join('\t') }
        }
        .flatten()
        .collectFile(name: 'name_replacement.txt', newLine: true)
        .ifEmpty([])
}

// Escape Python-regex metacharacters and YAML single-quote a sample ID
// for use in a multiqcSampleMergeYaml lookbehind pattern.
def multiqcSampleMergeYamlPattern(id, read) {
    def esc = id.replaceAll(/[\\^$.|?*+()\[\]{}\/]/) { m -> "\\${m[0]}" }
                .replace("'", "''")
    return "    - type: regex\n      pattern: '(?<=^${esc})_${read}\$'"
}

//
// MultiQC table_sample_merge YAML scoped to PE sample IDs via a
// fixed-length lookbehind, so sample IDs ending in `_1` / `_2` aren't
// wrongly collapsed.
//
def multiqcSampleMergeYaml(samplesheet_path, schema_path) {
    // Row order comes from assets/schema_input.json: [0]=meta,
    // [1]=fastq_1, [2]=fastq_2 (truthy => paired-end).
    def pe_sample_ids = samplesheetToList(samplesheet_path, schema_path)
        .findAll { row -> row[2] as boolean }
        .collect { row -> row[0].id as String }
        .unique()
        .sort()
    if (!pe_sample_ids) return 'table_sample_merge: {}\n'

    def r1 = pe_sample_ids.collect { id -> multiqcSampleMergeYamlPattern(id, 1) }.join('\n')
    def r2 = pe_sample_ids.collect { id -> multiqcSampleMergeYamlPattern(id, 2) }.join('\n')
    return "table_sample_merge:\n  \"Read 1\":\n${r1}\n  \"Read 2\":\n${r2}\n"
}

//
// Load a MultiQC custom-content config template from a YAML file. The
// asset is parsed as YAML so SnakeYAML stays contained in this single
// helper and callers just get a plain Map. Top-level keys starting
// with '_' are dropped after YAML anchor resolution (they exist only
// to host named anchors reused via merge keys elsewhere in the file),
// so they are never emitted to MultiQC.
//
def loadMultiqcAsset(asset_path) {
    def parsed = new org.yaml.snakeyaml.Yaml().load(file(asset_path).text)
    parsed.findAll { k, _v -> !k.toString().startsWith('_') }
}

//
// Certainty of a strand call, expressed as the same quantity
// `calculateStrandedness` compares against `stranded_threshold`: the
// inferred direction's share of the stranded fragment pool, 0-100.
// So a 'forward' call that cleared `stranded_threshold = 0.8` will
// show a value >= 80 here regardless of the unstranded fraction.
// Null when the input is null, when the sample has zero stranded
// fragments, or for 'unstranded' and 'undetermined' classifications
// (different thresholds apply to those calls).
//
def inferenceCertainty(analysis) {
    if (!analysis) return null
    def fwd = analysis.forwardFragments
    def rev = analysis.reverseFragments
    def stranded = fwd + rev
    if (stranded == 0) return null

    def s = analysis.inferred_strandedness
    if (s == 'forward') return (fwd / stranded) * 100
    if (s == 'reverse') return (rev / stranded) * 100
    null
}

// Round a Double to one decimal place, preserving null.
def roundOneDecimal(v) {
    v == null ? null : Math.round(v * 10) / 10.0d
}

//
// Build a per-sample cell map for the strandedness summary table. One
// entry per column id; nulls mean "method didn't produce this cell"
// and are dropped before emission so MultiQC renders blanks (not
// "None") in data exports.
//
def strandSummaryCells(_meta, provided, status, salmon, rseqc) {
    [
        provided:        provided,
        salmon_inferred: salmon?.inferred_strandedness ?: '-',
        salmon_pct:      roundOneDecimal(inferenceCertainty(salmon)),
        salmon_s:        roundOneDecimal(salmon?.forwardFragments),
        salmon_a:        roundOneDecimal(salmon?.reverseFragments),
        salmon_u:        roundOneDecimal(salmon?.unstrandedFragments),
        rseqc_inferred:  rseqc?.inferred_strandedness ?: '-',
        rseqc_pct:       roundOneDecimal(inferenceCertainty(rseqc)),
        rseqc_s:         roundOneDecimal(rseqc?.forwardFragments),
        rseqc_a:         roundOneDecimal(rseqc?.reverseFragments),
        rseqc_u:         roundOneDecimal(rseqc?.unstrandedFragments),
        status:          status,
    ]
}

//
// Build the MultiQC custom-content JSON for the strandedness summary
// table by merging a static config template (parsed from
// assets/strand_check_summary.yaml) with per-sample rows
// emitted by classifyStrand. Column order is taken from the YAML
// header keyset so reordering columns in the asset reorders them in
// the rendered table. Throws if a row emits a cell that is not
// declared in the asset's headers block, so the data/config contract
// stays explicit.
//
def strandCheckSummaryYaml(static_config, rows) {
    def header_keys = static_config.headers.keySet()
    // Sort by sample id so the merged output is deterministic regardless of
    // which sample finished RSeQC/Salmon first, and so the rendered MultiQC
    // table has a consistent default row order.
    def data = rows.toSorted { row -> row[0].id }.collectEntries { row ->
        def (meta, provided, status, salmon, rseqc) = row
        def raw = strandSummaryCells(meta, provided, status, salmon, rseqc)
        def unknown = raw.keySet() - header_keys
        if (unknown) error("strand_check_summary.yaml headers do not declare columns: ${unknown}")

        def cells = [:]  // follow header order, drop null cells
        header_keys.each { k -> if (raw[k] != null) cells[k] = raw[k] }
        [ (meta.id): cells ]
    }
    groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(static_config + [data: data]))
}

// Per-sample {Sense/Antisense/Unstranded} percentages for the strand
// composition bargraph. Returns null so unavailable datasets are
// dropped before rendering.
def strandCompositionMap(analysis) {
    if (!analysis) return null
    [
        Sense:      roundOneDecimal(analysis.forwardFragments),
        Antisense:  roundOneDecimal(analysis.reverseFragments),
        Unstranded: roundOneDecimal(analysis.unstrandedFragments),
    ]
}

//
// Build the MultiQC custom-content JSON for the strandedness read-
// composition bargraph. When both inference methods produced data,
// two datasets are emitted (RSeQC first so reports default to the
// alignment-based view) and MultiQC's `data_labels` switcher lets
// users flip between them. Single-dataset otherwise. Dataset labels
// inherit `ylab` from the static config's pconfig so the string lives
// in YAML only.
//
def strandCheckCompositionYaml(static_config, rows) {
    def rseqc_data  = [:]
    def salmon_data = [:]
    // Sort by sample id so the merged output is deterministic regardless of
    // which sample finished RSeQC/Salmon first, and so the rendered MultiQC
    // bargraph has a consistent default sample order.
    rows.toSorted { row -> row[0].id }.each { row ->
        def (meta, _p, _s, salmon, rseqc) = row
        if (rseqc)  rseqc_data[meta.id]  = strandCompositionMap(rseqc)
        if (salmon) salmon_data[meta.id] = strandCompositionMap(salmon)
    }
    def datasets = []
    def labels   = []
    if (rseqc_data)  { datasets << rseqc_data;  labels << 'RSeQC'  }
    if (salmon_data) { datasets << salmon_data; labels << 'Salmon' }

    // Deep-ish copy: both the top-level map and pconfig get mutated,
    // so clone both to keep the cached static_config untouched.
    def config  = new LinkedHashMap(static_config)
    def pconfig = new LinkedHashMap(config.pconfig)
    if (datasets.size() > 1) {
        pconfig.data_labels = labels.collect { label -> [name: label, ylab: pconfig.ylab] }
    }
    config.pconfig = pconfig
    config.data    = datasets.size() == 1 ? datasets[0] : datasets
    groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
}
