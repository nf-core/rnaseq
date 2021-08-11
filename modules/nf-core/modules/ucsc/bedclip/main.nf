// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '377'

process UCSC_BEDCLIP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ucsc-bedclip=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bedclip:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2"
    }

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedClip \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        $software: \$(echo $VERSION)
    END_VERSIONS
    """
}
