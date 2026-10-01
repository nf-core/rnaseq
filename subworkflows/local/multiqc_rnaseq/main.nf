nextflow.enable.types = true

//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                    } from '../../../modules/nf-core/multiqc'
include { workflowVersionToYAML      } from '../../nf-core/utils_nfcore_pipeline'
include { Sample; MultiqcFiles; SampleRuns; TrimReadCount; PercentMappedPass; StrandData } from '../../../modules/nf-core/types'
include { MultiqcReport } from '../../../modules/nf-core/multiqc/main'
include { methodsDescriptionText     } from '../utils_nfcore_rnaseq_pipeline'
include { workflowSummaryMultiqcYaml } from './helpers'
include { multiqcNameReplacementLines } from './helpers'
include { multiqcSampleMergeYaml     } from './helpers'
include { loadMultiqcAsset           } from './helpers'
include { strandCheckSummaryYaml     } from './helpers'
include { strandCheckCompositionYaml } from './helpers'

workflow MULTIQC_RNASEQ {

    take:
    ch_sample_ids: Channel<String>                // one id per input sample; every sample gets a report under skip_quantification_merge
    ch_mqc_files: Channel<MultiqcFiles>           // per-sample files from each stage, for both report modes
    ch_mqc_sample_only: Channel<MultiqcFiles>     // per-sample files for the per-sample reports only
    ch_mqc_report_only: Channel<Path>             // files for the merged report only
    ch_strand_data: Channel<StrandData>           // per-sample strand classification, used for the Strandedness checks section
    ch_trim_read_count: Channel<TrimReadCount>  // for fail_trimmed section
    ch_percent_mapped_pass: Channel<PercentMappedPass> // for fail_mapped section
    aligner_display_name: String                  // display name of the aligner used for the percent_mapped metric, e.g. 'STAR uniquely mapped reads' or 'Bowtie2 overall alignment rate'
    ch_fastq: Channel<SampleRuns>                 // one entry per sample, one run per sequencing run
    ch_collated_versions: Channel<Path>           // versions yaml
    mqc_default_config: Path                      // pipeline-bundled MultiQC config
    mqc_custom_config: Path?                      // optional user MultiQC config
    mqc_logo: Path?                               // optional custom logo
    methods_description_yml: Path                 // methods-description YAML template
    strand_summary_asset: Path                    // strand_check_summary YAML custom-content template
    strand_composition_asset: Path                // strand_check_composition YAML custom-content template
    sample_status_header: Path                    // MultiQC custom content header for fail_* tables
    min_trimmed_reads: Integer                    // threshold for fail_trimmed classification
    skip_quantification_merge: Boolean

    main:

    //
    // fail_* custom-content TSVs. Each failing sample contributes one
    // file; the merged tables concatenate them in sample order and are
    // only produced when at least one sample fails.
    //
    // `status_header_lines` tracks the header row count so editing
    // `sample_status_header.txt` doesn't silently mis-skip the merged
    // aggregate's concatenation.
    //
    status_header_lines = sample_status_header.readLines().size() + 1  // parent header + one column row
    status_header_text  = sample_status_header.text

    ch_fail_trimmed_rows = ch_trim_read_count
        .filter { r -> r.num_reads <= min_trimmed_reads }
        .map { r ->
            record(
                id:      r.id,
                name:    "${r.id}_fail_trimmed_samples_mqc.tsv",
                content: "Sample\tReads after trimming\n${r.id}\t${r.num_reads}\n"
            )
        }

    ch_fail_trimmed_by_id = ch_fail_trimmed_rows
        .collectFile { r -> [r.name, r.content] }
        .map { p -> p as Path }
        .map { f -> record(id: f.name.replace('_fail_trimmed_samples_mqc.tsv', ''), fail_trimmed: f) }

    ch_fail_trimmed_merged = ch_fail_trimmed_rows
        .map { r -> r.content }
        .collectFile(name: 'fail_trimmed_samples_mqc.tsv', keepHeader: true, sort: true)
        .map { p -> p as Path }

    ch_fail_mapped_rows = ch_percent_mapped_pass
        .filter { r -> r.pass != null && !r.pass }
        .map { r ->
            record(
                id:      r.id,
                name:    "${r.id}_fail_mapped_samples_mqc.tsv",
                content: status_header_text + "Sample\t${aligner_display_name} (%)\n${r.id}\t${r.percent_mapped}\n"
            )
        }

    ch_fail_mapped_by_id = ch_fail_mapped_rows
        .collectFile { r -> [r.name, r.content] }
        .map { p -> p as Path }
        .map { f -> record(id: f.name.replace('_fail_mapped_samples_mqc.tsv', ''), fail_mapped: f) }

    ch_fail_mapped_merged = ch_fail_mapped_rows
        .map { r -> r.content }
        .collectFile(name: 'fail_mapped_samples_mqc.tsv', keepHeader: true, skip: status_header_lines, sort: true)
        .map { p -> p as Path }

    //
    // Strandedness checks custom-content section. Two MultiQC
    // subsections (summary table + stacked composition bargraph) are
    // rendered from the same per-sample record, with header / pconfig
    // / colour config in the bundled YAML templates. The composition
    // section inherits `parent_*` from the summary section so the
    // description lives in one place.
    //
    strand_summary_static     = loadMultiqcAsset(strand_summary_asset)
    strand_composition_static = loadMultiqcAsset(strand_composition_asset) + strand_summary_static.subMap(['parent_id', 'parent_name', 'parent_description'])

    // Per-run table_sample_merge config: only PE samples from the
    // input get their _1 / _2 rows grouped in the General Stats
    // table.
    ch_mqc_dynamic_config = ch_fastq
        .collect()
        .flatMap { samples -> [ multiqcSampleMergeYaml(samples) ] }
        .collectFile(name: 'multiqc_sample_merge.yml')
        .collect()
        .map { files -> files.toList().first() as Path }

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel.of(workflowSummaryMultiqcYaml())
        .collectFile(name: 'workflow_summary_mqc.yaml')
        .collect()
        .map { files -> files.toList().first() as Path }

    ch_methods_description = channel.of(methodsDescriptionText(methods_description_yml))
        .collectFile(name: 'methods_description_mqc.yaml')
        .collect()
        .map { files -> files.toList().first() as Path }

    //
    // Two execution modes for MULTIQC:
    //   - merged (default): one report covers the whole run.
    //   - per-sample (--skip_quantification_merge): one report per
    //     sample; workflow-level versions are replaced with a
    //     pipeline-identity manifest so the report doesn't wait on
    //     the global versions topic.
    //
    // Each branch ends with a record matching the MULTIQC input
    // contract (id, meta, files, configs, logo, replace_names, sample_names).
    //
    if (skip_quantification_merge) {
        ch_strand_summary_by_id = ch_strand_data
            .collectFile { r -> ["${r.id}_strand_check_summary_mqc.json", strandCheckSummaryYaml(strand_summary_static, [r]) as String] }
            .map { p -> p as Path }
            .map { f -> record(id: f.name.replace('_strand_check_summary_mqc.json', ''), strand_summary: f) }

        ch_strand_composition_by_id = ch_strand_data
            .collectFile { r -> ["${r.id}_strand_check_composition_mqc.json", strandCheckCompositionYaml(strand_composition_static, [r]) as String] }
            .map { p -> p as Path }
            .map { f -> record(id: f.name.replace('_strand_check_composition_mqc.json', ''), strand_composition: f) }

        // One empty contribution per sample keeps samples that no stage contributed files for.
        def ch_no_files: Channel<MultiqcFiles> = ch_sample_ids.map { id -> record(id: id, files: []) }
        ch_per_sample_bundle = ch_no_files
            .mix(ch_mqc_files)
            .mix(ch_mqc_sample_only)
            .mix(ch_fail_trimmed_by_id.map { r -> record(id: r.id, files: [r.fail_trimmed]) })
            .mix(ch_fail_mapped_by_id.map { r -> record(id: r.id, files: [r.fail_mapped]) })
            .mix(ch_strand_summary_by_id.map { r -> record(id: r.id, files: [r.strand_summary]) })
            .mix(ch_strand_composition_by_id.map { r -> record(id: r.id, files: [r.strand_composition]) })
            .collect()
            .flatMap { rs ->
                rs.collect { r -> r.id }.toSet().toSorted().collect { id ->
                    def files = rs.findAll { r -> r.id == id }.collectMany { r -> r.files }
                    record(id: id, files: files.toSorted { f -> f.name })
                }
            }

        ch_manifest_versions = channel.of(workflowVersionToYAML())
            .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml')
            .collect()
            .map { files -> files.toList().first() as Path }

        ch_static_globals = ch_workflow_summary
            .combine(ch_methods_description)
            .combine(ch_manifest_versions)
            .map { workflow_summary, methods_description, manifest_versions -> [workflow_summary, methods_description, manifest_versions] }

        ch_global_files = ch_fail_trimmed_merged
            .mix(ch_fail_mapped_merged)
            .collect()
            .map { files -> files.toSorted { f -> f.name } }

        ch_multiqc_input = ch_per_sample_bundle
            .combine(static_globals: ch_static_globals, run_globals: ch_global_files, dyn: ch_mqc_dynamic_config)
            .map { r ->
                // No replace_names: each per-sample report contains one sample.
                record(
                    id:             r.id,
                    meta:           [id: r.id],
                    multiqc_files:  r.files + r.static_globals + r.run_globals,
                    multiqc_config: [mqc_default_config, r.dyn, mqc_custom_config].findAll { cfg -> cfg != null }.toList(),
                    multiqc_logo:   mqc_logo,
                    replace_names:  null,
                    sample_names:   null
                )
            }
    } else {
        // Zero strand rows -> no *_mqc.json emission -> MultiQC drops
        // the section cleanly.
        ch_strand_rows = ch_strand_data.collect()

        ch_strand_summary_merged = ch_strand_rows
            .flatMap { rows -> rows.isEmpty() ? [] : [strandCheckSummaryYaml(strand_summary_static, rows)] }
            .collectFile(name: 'strand_check_summary_mqc.json')
            .map { p -> p as Path }

        ch_strand_composition_merged = ch_strand_rows
            .flatMap { rows -> rows.isEmpty() ? [] : [strandCheckCompositionYaml(strand_composition_static, rows)] }
            .collectFile(name: 'strand_check_composition_mqc.json')
            .map { p -> p as Path }

        // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
        ch_name_replacements = ch_fastq
            .collect()
            .flatMap { rows -> multiqcNameReplacementLines(rows) }
            .collectFile(name: 'name_replacement.txt', newLine: true)
            .map { p -> p as Path }
            .collect()
            .map { files -> files.isEmpty() ? null : files.toList().first() }

        // `multiqc_report` is a sentinel meta.id used by
        // conf/modules/multiqc.config to pick the merged output path.
        ch_multiqc_files_merged = ch_mqc_files
            .flatMap { r -> r.files }
            .mix(ch_mqc_report_only)
            .mix(ch_fail_trimmed_merged)
            .mix(ch_fail_mapped_merged)
            .mix(ch_strand_summary_merged)
            .mix(ch_strand_composition_merged)
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        ch_multiqc_input = ch_multiqc_files_merged
            .collect()
            .map { files -> record(files: files.toSorted { f -> f.name }) }
            .combine(replace_names: ch_name_replacements, dyn: ch_mqc_dynamic_config)
            .flatMap { r ->
                [
                    record(
                        id:             'multiqc_report',
                        meta:           [id: 'multiqc_report'],
                        multiqc_files:  r.files,
                        multiqc_config: [mqc_default_config, r.dyn, mqc_custom_config].findAll { cfg -> cfg != null }.toList(),
                        multiqc_logo:   mqc_logo,
                        replace_names:  r.replace_names,
                        sample_names:   null
                    )
                ]
            }
    }

    //
    // One record per MULTIQC task: a single 'multiqc_report' row when
    // merged, or one per sample under skip_quantification_merge.
    //
    ch_results = MULTIQC(ch_multiqc_input)

    emit:
    ch_results // channel: MultiqcReport
}
