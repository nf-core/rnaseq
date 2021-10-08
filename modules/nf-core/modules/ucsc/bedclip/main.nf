// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

def VERSION = '377'

process UCSC_BEDCLIP {
    tag "$meta.id"
    label 'process_medium'

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
    path "versions.yml"                , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    bedClip \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
