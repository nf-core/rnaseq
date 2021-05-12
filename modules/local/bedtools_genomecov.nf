// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sense.bedGraph")    , emit: bedgraph_sense
    tuple val(meta), path("*.antisense.bedGraph"), emit: bedgraph_antisense
    path "*.version.txt"                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def prefix_sense     = "${prefix}.sense"
    def prefix_antisense = "${prefix}.antisense"
    if (meta.strandedness == 'reverse') {
        prefix_sense     = "${prefix}.antisense"
        prefix_antisense = "${prefix}.sense"
    }
    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand + \\
        $options.args \\
        | bedtools sort > ${prefix_sense}.bedGraph

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand - \\
        $options.args \\
        | bedtools sort > ${prefix_antisense}.bedGraph

    bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
    """
}
