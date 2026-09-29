nextflow.enable.types = true

process FQ_LINT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'quay.io/biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    tuple(meta: Map, fastq: List<Path>)

    output:
    tuple(meta, file("*.fq_lint.txt"))

    topic:
    tuple(task.process, 'fq', eval("fq lint --version | sed 's/fq-lint //; s/ .*//'")) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq lint \\
        $args \\
        ${fastq.join(' ')} > ${prefix}.fq_lint.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fq_lint.txt
    """
}
