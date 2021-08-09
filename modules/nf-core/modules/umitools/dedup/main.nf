// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0"
    } else {
        container "quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path   "versions.yml"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def paired   = meta.single_end ? "" : "--paired"
    """
    umi_tools dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        $paired \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(umi_tools --version 2>&1) | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
