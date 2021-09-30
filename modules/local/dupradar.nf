// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DUPRADAR {
    tag "$meta.id"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-dupradar=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.18.0--r40_1"
    } else {
        container "quay.io/biocontainers/bioconductor-dupradar:1.18.0--r40_1"
    }

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.txt")    , emit: txt
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path "versions.yml"               , emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired_end = meta.single_end ? 'single' :  'paired'
    """
    dupradar.r \\
        $bam \\
        $prefix \\
        $gtf \\
        $strandedness \\
        $paired_end \\
        $task.cpus

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dupradar: \$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
    END_VERSIONS
    """
}
