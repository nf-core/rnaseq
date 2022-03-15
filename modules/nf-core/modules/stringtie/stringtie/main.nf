process STRINGTIE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0' :
        'quay.io/biocontainers/stringtie:2.1.7--h978d192_0' }"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.coverage.gtf")   , emit: coverage_gtf
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), path("*.abundance.txt")  , emit: abundance
    tuple val(meta), path("*.ballgown")       , emit: ballgown
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

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
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
