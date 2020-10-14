// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SAMTOOLS_INDEX {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bai"), emit: bai
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    samtools index $bam
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
