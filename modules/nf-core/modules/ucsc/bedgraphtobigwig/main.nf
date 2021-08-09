// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '377'

process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

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
    path  "versions.yml"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedGraphToBigWig $bedgraph $sizes ${prefix}.bigWig
    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - $software: \$(echo $VERSION)
    END_VERSIONS
    """
}
