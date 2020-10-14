// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSION = '2.2.0'

process HISAT2_EXTRACTSPLICESITES {
    tag "$gtf"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::hisat2=2.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hisat2:2.2.0--py37hfa133b6_4"
    } else {
        container "quay.io/biocontainers/hisat2:2.2.0--py37hfa133b6_4"
    }

    input:
    path gtf

    output:
    path "*.splice_sites.txt", emit: txt
    path "*.version.txt"     , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.splice_sites.txt
    echo $VERSION > ${software}.version.txt
    """
}
