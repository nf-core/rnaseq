// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def VERSION = '377'

process UCSC_BEDRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    //container "https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1"

    conda (params.conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)

    input:
    tuple val(meta), path(bedgraph)
    path sizes
    val options

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    bedGraphToBigWig $bedgraph $sizes ${prefix}.bigWig
    echo $VERSION > ${software}.version.txt
    """
}
