// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
options              = initOptions(params.options)

def VERSION = '1.38.0'

process DEXSEQ_ANNOTATE {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::htseq=0.13.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/htseq:0.13.5--py39h70b41aa_1"
    } else {
        container "quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1"
    }

    input:
    path gtf

    output:
    path "dexseq.gff"             , emit: dexseq_gff
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    dexseq_prepare_annotation.py $gtf dexseq.gff
    echo $VERSION > ${software}.version.txt
    """
}
