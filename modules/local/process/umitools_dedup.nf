// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/umi_tools:1.0.1--py37h516909a_1"
    //container "https://depot.galaxyproject.org/singularity/umi_tools:1.0.1--py37h516909a_1"

    conda (params.conda ? "bioconda::umi_tools=1.0.1" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    umi_tools dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        --output-stats=$prefix \\
        $ioptions.args \\

    umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
    """
}
