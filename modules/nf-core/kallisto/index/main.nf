process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("kallisto")  , emit: index
    tuple val("${task.process}"), val('kallisto'), cmd("echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//'"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto \\
        index \\
        $args \\
        -i kallisto \\
        $fasta
    """

    stub:
    """
    touch kallisto
    """
}
