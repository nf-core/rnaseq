nextflow.enable.types = true

record MultiqcInput {
    id:             String
    meta:           Map
    multiqc_files:  List<Path>
    multiqc_config: List<Path>
    multiqc_logo:   Path?
    replace_names:  Path?
    sample_names:   Path?
}

process MULTIQC {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    sample: MultiqcInput

    stage:
    stageAs sample.multiqc_files, '?/*'
    stageAs sample.multiqc_config, '?/*'

    // MultiQC must not push its version to the `versions` topic: its input depends on that topic, so the pipeline would hang forever
    output:
    record(
        id:     sample.id,
        meta:   sample.meta,
        report: file('*.html'),
        data:   file('*_data'),
        plots:  file('*_plots', optional: true)
    )

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = sample.multiqc_config ? sample.multiqc_config instanceof List ? "--config ${sample.multiqc_config.join(' --config ')}" : "--config ${sample.multiqc_config}" : ""
    def logo = sample.multiqc_logo ? "--cl-config 'custom_logo: \"${sample.multiqc_logo}\"'" : ''
    def replace = sample.replace_names ? "--replace-names ${sample.replace_names}" : ''
    def samples = sample.sample_names ? "--sample-names ${sample.sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}
