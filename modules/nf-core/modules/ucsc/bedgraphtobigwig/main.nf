def VERSION = '377'

process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1"
    } else {
        container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    }

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    path "versions.yml"              , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bigWig

    cat <<-END_VERSIONS > versions.yml
    UCSC_BEDGRAPHTOBIGWIG:
        ucsc: \$(echo $VERSION)
    END_VERSIONS
    """
}
