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
    tuple val(meta), path("*.forward.bedGraph"), emit: bedgraph_forward
    tuple val(meta), path("*.reverse.bedGraph"), emit: bedgraph_reverse
    path  "versions.yml"                       , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def prefix_forward = "${prefix}.forward"
    def prefix_reverse = "${prefix}.reverse"
    if (meta.strandedness == 'reverse') {
        prefix_forward = "${prefix}.reverse"
        prefix_reverse = "${prefix}.forward"
    }
    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand + \\
        $options.args \\
        | bedtools sort > ${prefix_forward}.bedGraph

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand - \\
        $options.args \\
        | bedtools sort > ${prefix_reverse}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
