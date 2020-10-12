// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEQC_BAMSTAT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0"
    //container "https://depot.galaxyproject.org/singularity/rseqc:4.0.0--py38h0213d0e_0"

    conda (params.conda ? "bioconda::rseqc=4.0.0" : null)

    input:
    tuple val(meta), path(bam)
    val   options

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt
    path  "*.version.txt"                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    bam_stat.py \\
        -i $bam \\
        $ioptions.args \\
        > ${prefix}.bam_stat.txt

    bam_stat.py --version | sed -e "s/bam_stat.py //g" > ${software}.version.txt
    """
}
