// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SALMON_QUANT {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::salmon=1.3.0" : null)
    if (workflow.containerEngine == 'singularity' || !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.3.0--hf69c8f4_0"
    } else {
        container "quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    
    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software  = getSoftwareName(task.process)
    prefix        = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def endedness = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"

    def strandedness = meta.single_end ? 'U' : 'IU'
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? 'SF' : 'ISF'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? 'SR' : 'ISR'
    }
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        --index $index \\
        $endedness \\
        $options.args \\
        -o $prefix

    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}
