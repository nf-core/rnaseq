// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '1.38.0'

process DEXSEQ_COUNT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::htseq=0.13.5 bioconda::pysam=0.16.0.1-3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-99b86c9c0de90c62fa0fb4ad6daef5bf8b261885:8e553aed88d367cdb3ea41614ad895bae8c9818b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-99b86c9c0de90c62fa0fb4ad6daef5bf8b261885:8e553aed88d367cdb3ea41614ad895bae8c9818b-0"
    }

    input:
    tuple val(meta), path(bam)
    path dexseq_gff

    output:
    tuple val(meta), path("*.counts.txt"), emit: dexseq_counts
    path "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    dexseq_count.py \\
        -f bam \\
        $options.args \\
        $dexseq_gff \\
        $bam \\
        ${prefix}.txt

    echo $VERSION > ${software}.version.txt
    """
}