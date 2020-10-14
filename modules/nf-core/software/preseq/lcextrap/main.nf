// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::preseq=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/preseq:2.0.3--hf53bd2b_3"
    } else {
        container "quay.io/biocontainers/preseq:2.0.3--hf53bd2b_3"
    }

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log")       , emit: log
    path  "*.version.txt"                , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $options.args \\
        $paired_end \\
        -output ${prefix}.ccurve.txt \\
        $bam
    cp .command.err ${prefix}.command.log

    echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//' > ${software}.version.txt
    """
}
