// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::umi_tools=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/umi_tools:1.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/umi_tools:1.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.tsv"), emit: tsv
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    umi_tools dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        --output-stats=$prefix \\
        $options.args \\

    umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
    """
}
