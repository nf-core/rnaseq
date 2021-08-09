// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::preseq=3.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/preseq:3.1.2--h06ef8b0_1"
    } else {
        container "quay.io/biocontainers/preseq:3.1.2--h06ef8b0_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log")       , emit: log
    path   "versions.yml"                , emit: version

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

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
