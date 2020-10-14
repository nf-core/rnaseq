// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    //container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"

    conda (params.conda ? "bioconda::rseqc=3.0.1" : null)

    input:
    tuple val(meta), path(bam)
    path  bed
    
    output:
    tuple val(meta), path("*.junction.bed"), emit: bed
    tuple val(meta), path("*.Interact.bed"), emit: interact_bed
    tuple val(meta), path("*.xls")         , emit: xls
    tuple val(meta), path("*junction.pdf") , emit: pdf
    tuple val(meta), path("*events.pdf")   , emit: events_pdf
    tuple val(meta), path("*.r")           , emit: rscript
    tuple val(meta), path("*.log")         , emit: log
    path  "*.version.txt"                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $options.args \\
        2> ${prefix}.junction_annotation.log

    junction_annotation.py --version | sed -e "s/junction_annotation.py //g" > ${software}.version.txt
    """
}
