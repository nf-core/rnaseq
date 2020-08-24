// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    //container " https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val options

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    samtools flagstat $bam > ${bam}.flagstat
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
