// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUNZIP {
    tag "$archive"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "biocontainers/biocontainers:v1.2.0_cv1"

    conda (params.conda ? "conda-forge::sed=4.7" : null)
    
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
