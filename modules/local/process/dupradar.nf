// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DUPRADAR {
    tag "$meta.id"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam)
    path  gtf
    val   options

    output:
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.txt")    , emit: txt
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path  "*.version.txt"             , emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired_end = meta.single_end ? 'single' :  'paired'
    """
    dupradar.r $bam $gtf $strandedness $paired_end $task.cpus
    Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='${software}.version.txt')"
    """
}
