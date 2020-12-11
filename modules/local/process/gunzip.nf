// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUNZIP {
    tag "$archive"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    // conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    // } else {
    //     container "biocontainers/biocontainers:v1.2.0_cv1"
    // }
    
    input:
    path archive
    
    output:
    path "$gunzip",       emit: gunzip
    path "*.version.txt", emit: version
    
    script:
    def software = getSoftwareName(task.process)
    gunzip       = archive.toString() - '.gz'
    """
    gunzip --force $options.args $archive
    echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//' > ${software}.version.txt
    """
}
