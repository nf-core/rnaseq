// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
options              = initOptions(params.options)

def VERSION = '1.38.0'

process DEXSEQ_ANNOTATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-dexseq=1.38.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dexseq%3A1.38.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dexseq:1.38.0--r41hdfd78af_0"
    }

    input:
    path gtf

    output:
    path "dexseq.gff"             , emit: dexseq_gff
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    dexseq_prepare_annotation.py ${options.args} $gtf dexseq.gff
    echo $VERSION > ${software}.version.txt
    """
}
