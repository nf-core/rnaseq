// Import generic module functions
include { saveFiles; getModuleName } from './functions'

params.options = [:]

process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    cache false

    input:
    path versions

    output:
    path "software_versions.yml"     , emit: yml
    path 'software_versions_mqc.yaml', emit: mqc_yaml

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    """
    cat - $versions <<-END_WORKFLOW_VERSION > software_versions.yml
    Workflow:
        - Nextflow: $workflow.nextflow.version
        - $workflow.manifest.name: $workflow.manifest.version
    END_WORKFLOW_VERSION
    cat - <<-END_MQC_YAML > software_versions_mqc.yaml
    id: 'software_versions'
    section_name: '{{ name }} Software Versions'
    section_href: 'https://github.com/{{ name }}'
    plot_type: 'table'
    description: 'are collected at run time from the software output.'
    data:
    \$( sed 's/^/    /' software_versions.yml )
    END_MQC_YAML
    """
}
