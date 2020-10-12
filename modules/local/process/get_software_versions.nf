// Import generic module functions
include { saveFiles } from './functions'

/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'pipeline_info', publish_id:'') }

    container "quay.io/biocontainers/python:3.8.3"
    //container  https://depot.galaxyproject.org/singularity/python:3.8.3

    conda (params.conda ? "conda-forge::python=3.8.3" : null)

    input:
    path versions
    val  options

    output:
    path "software_versions.csv"     , emit: csv
    path 'software_versions_mqc.yaml', emit: yaml

    script:
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
