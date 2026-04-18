//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                 } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcNameReplacements } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcSampleMergeYaml  } from '../utils_nfcore_rnaseq_pipeline'

workflow MULTIQC_RNASEQ {

    take:
    ch_multiqc_files                  // channel: [ val(meta), path(file) ]    - flat, merged-mode input
    ch_per_sample_bundle              // channel: [ val(meta), [ files ] ]     - per-sample mode, one tuple per sample
    ch_per_sample_globals             // channel: [ [:], path(file) ]          - per-sample mode, run-level files
    ch_per_sample_collated_versions   // channel: path(versions yaml)          - per-sample mode, early-closing
    ch_fastq                          // channel: [ val(meta), [ reads ] ]
    ch_collated_versions              // channel: path(versions yaml)          - merged mode
    samplesheet_path                  // path: pipeline input samplesheet
    samplesheet_schema                // path: samplesheet JSON schema
    mqc_default_config                // path: pipeline-bundled MultiQC config
    mqc_custom_config                 // path (or []): optional user MultiQC config
    mqc_logo                          // path (or []): optional custom logo
    methods_description_yml           // path: methods-description YAML template
    skip_quantification_merge         // boolean

    main:

    // Per-run table_sample_merge config: only PE samples from the samplesheet
    // get their _1/_2 rows grouped in the General Stats table.
    ch_mqc_dynamic_config = channel.of(multiqcSampleMergeYaml(samplesheet_path, samplesheet_schema))
        .collectFile(name: 'multiqc_sample_merge.yml')

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel.value(
        paramsSummaryMultiqc(
            paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        )
    ).collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description = channel.value(
        methodsDescriptionText(methods_description_yml)
    ).collectFile(name: 'methods_description_mqc.yaml')

    // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
    ch_name_replacements = multiqcNameReplacements(ch_fastq)

    if (skip_quantification_merge) {
        // One MultiQC report per sample. Each incoming tuple is already bundled
        // as [meta, files_list] for a single sample — no groupTuple, no count
        // helper. Globals (workflow summary, methods description, run-level
        // stragglers like fail_trimmed_samples_mqc, and a per-sample-closing
        // versions YAML) are broadcast to every per-sample report via `combine`.
        ch_static_globals = ch_workflow_summary
            .mix(ch_methods_description)
            .mix(ch_per_sample_collated_versions)
            .collect()

        ch_run_level_globals = ch_per_sample_globals
            .map { _meta, f -> f }
            .collect()
            .ifEmpty([])

        ch_multiqc_input = ch_per_sample_bundle
            .combine(ch_static_globals.toList())
            .combine(ch_run_level_globals.toList())
            .combine(ch_mqc_dynamic_config)
            .map { meta, sample_files, static_globals, run_globals, dyn ->
                def all_globals = (static_globals ?: []) + (run_globals ?: [])
                [
                    [id: meta.id],
                    sample_files + all_globals,
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
                    mqc_logo,
                    [],  // no replace_names - each report contains one sample's files
                    [],
                ]
            }
    } else {
        // One merged MultiQC report. 'multiqc_report' is a sentinel meta.id
        // used by conf/modules/multiqc.config to pick the merged output
        // path/prefix. Wrap the collected file list in a 1-tuple so
        // .combine() doesn't spread it across the downstream closure args.
        ch_multiqc_all = ch_multiqc_files.mix(
            ch_workflow_summary
                .mix(ch_collated_versions)
                .mix(ch_methods_description)
                .map { f -> [[:], f] }
        )

        ch_all_files = ch_multiqc_all
            .map { _meta, f -> f }
            .collect()
            .map { files -> [files] }

        ch_multiqc_input = ch_all_files
            .combine(ch_name_replacements.ifEmpty([]).toList())
            .combine(ch_mqc_dynamic_config)
            .map { files, replace_names, dyn ->
                [
                    [id: 'multiqc_report'],
                    files,
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
                    mqc_logo,
                    replace_names ?: [],
                    [],
                ]
            }
    }

    MULTIQC(ch_multiqc_input)

    emit:
    report = MULTIQC.out.report.map { _meta, report -> report }
}
