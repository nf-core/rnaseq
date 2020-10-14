// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SAMTOOLS_STATS {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity') {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "*.version.txt"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    samtools stats $bam > ${bam}.stats
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
