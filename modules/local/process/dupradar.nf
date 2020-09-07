// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DUPRADAR {
    tag "$meta.id"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // container "quay.io/biocontainers/bioconductor-dupradar:1.18.0--r40_0"
    // //container "https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.18.0--r40_0"
    //
    // conda (params.conda ? "bioconda::bioconductor-dupradar=1.18.0" : null)
    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(bam)
    path  gtf
    val   options

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.txt"), emit: txt
    path  "*.version.txt"         , emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    // Figure out strandedness from pipeline parameters
    def unstranded       = params.unstranded
    def forward_stranded = params.forward_stranded
    def reverse_stranded = params.reverse_stranded
    if (params.pico) {
        unstranded       = false
        forward_stranded = true
        reverse_stranded = false
    }
    def strandedness = 0
    if (forward_stranded && !unstranded) {
        strandedness = 1
    } else if (reverse_stranded && !unstranded) {
        strandedness = 2
    }
    def paired_end = meta.single_end ? 'single' :  'paired'
    """
    dupradar.r $bam $gtf $strandedness $paired_end $task.cpus
    Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='${software}.version.txt')"
    """
}
