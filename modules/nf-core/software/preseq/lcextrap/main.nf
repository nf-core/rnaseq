// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/preseq:2.0.3--hf53bd2b_3"
    //container "https://depot.galaxyproject.org/singularity/preseq:2.0.3--hf53bd2b_3"

    conda (params.conda ? "bioconda::preseq=2.0.3" : null)

    input:
    tuple val(meta), path(bam)
    val options

    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def pe       = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $ioptions.args \\
        $pe \\
        -output ${prefix}.ccurve.txt \\
        $bam
    cp .command.err ${prefix}.command.log

    echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//' > ${software}.version.txt
    """
}
