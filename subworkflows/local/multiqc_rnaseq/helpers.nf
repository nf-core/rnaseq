//
// Untyped helper functions for MULTIQC_RNASEQ. They live outside the typed
// subworkflow script because `nextflow lint` rejects calls to nf-schema plugin functions in typed
// scripts (nextflow-io/nextflow#7720) and because they use dynamic maps.
//

include { paramsSummaryMap     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Workflow summary rendered as a MultiQC custom-content YAML section.
//
def workflowSummaryMultiqcYaml() {
    return paramsSummaryMultiqc(paramsSummaryMap(workflow, parameters_schema: 'nextflow_schema.json'))
}

//
// MultiQC `--replace-names` lines: map each FASTQ simpleName to
// '<id>_1' / '<id>_2' (or '<id>' for SE), skipping cases where the
// simpleName already equals the sample ID (see #1341 / #1659).
// `fastq_rows` is a list of record(id, meta, runs), runs being [ [ fastq_1, fastq_2? ], ... ].
//
def multiqcNameReplacementLines(fastq_rows) {
    def lines = fastq_rows.collectMany { row ->
        def meta  = row.meta
        def reads = row.runs
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
    return lines.sort()
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
def multiqcSampleMergeYaml(samplesheet_rows) {
    // Row order comes from assets/schema_input.json: [0]=meta,
    // [1]=fastq_1, [2]=fastq_2 (truthy => paired-end).
    def pe_sample_ids = samplesheet_rows
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
    def data = rows.toSorted { row -> row.id }.collectEntries { row ->
        def raw = strandSummaryCells(row.meta, row.provided, row.status, row.salmon, row.rseqc)
        def unknown = raw.keySet() - header_keys
        if (unknown) error("strand_check_summary.yaml headers do not declare columns: ${unknown}")

        def cells = [:]  // follow header order, drop null cells
        header_keys.each { k -> if (raw[k] != null) cells[k] = raw[k] }
        [ (row.id): cells ]
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
    rows.toSorted { row -> row.id }.each { row ->
        if (row.rseqc)  rseqc_data[row.id]  = strandCompositionMap(row.rseqc)
        if (row.salmon) salmon_data[row.id] = strandCompositionMap(row.salmon)
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
