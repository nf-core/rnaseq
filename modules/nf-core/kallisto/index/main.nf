process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--heb0cbe2_0':
        'biocontainers/kallisto:0.51.1--heb0cbe2_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("kallisto")  , emit: index
    tuple val("${task.process}"), val('kallisto'), eval('kallisto 2>&1 | head -1 | sed "s/^kallisto //; s/Usage.*//"'), emit: versions_kallisto, topic: versions

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
    mkdir kallisto
    """
}
