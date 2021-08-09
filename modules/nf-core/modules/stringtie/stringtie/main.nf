// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRINGTIE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"
    }

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.coverage.gtf")   , emit: coverage_gtf
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), path("*.abundance.txt")  , emit: abundance
    tuple val(meta), path("*.ballgown")       , emit: ballgown
    path   "versions.yml"                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandedness \\
        -G $gtf \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene.abundance.txt \\
        -C ${prefix}.coverage.gtf \\
        -b ${prefix}.ballgown \\
        -p $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(stringtie --version 2>&1))
    END_VERSIONS
    """
}
