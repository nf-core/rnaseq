// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GFFREAD {
    tag "$gff"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0"
    } else {
        container "quay.io/biocontainers/gffread:0.12.1--h8b12597_0"
    }

    input:
    path gff

    output:
    path "*.gtf"        , emit: gtf
    path "versions.yml" , emit: versions

    script:
    def prefix   = options.suffix ? "${gff.baseName}${options.suffix}" : "${gff.baseName}"
    """
    gffread \\
        $gff \\
        $options.args \\
        -o ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
