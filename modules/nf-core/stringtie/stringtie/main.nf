process STRINGTIE_STRINGTIE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0' :
        'biocontainers/stringtie:2.2.3--h43eeafb_0' }"

    input:
    tuple val(meta), path(bam)
    path  annotation_gtf

    output:
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), path("*.abundance.txt")  , emit: abundance
    tuple val(meta), path("*.coverage.gtf")   , optional: true, emit: coverage_gtf
    tuple val(meta), path("*.ballgown")       , optional: true, emit: ballgown
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    def ballgown  = annotation_gtf ? "-b ${prefix}.ballgown" : ""
    def coverage  = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""

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
        $reference \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene.abundance.txt \\
        $coverage \\
        $ballgown \\
        -p $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.gene.abundance.txt
    touch ${prefix}.coverage.gtf
    touch ${prefix}.ballgown

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
