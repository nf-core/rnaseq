nextflow.enable.types = true


record FqLintInput {
    id:    String
    meta:  Map
    reads: List<Path>
    args:  String?
}

record FqLintResult {
    id:   String
    meta: Map
    lint: Path
}

process FQ_LINT {
    tag "$sample.meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'quay.io/biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    sample: FqLintInput

    output:
    record(id: sample.id, meta: sample.meta, lint: file("*.fq_lint.txt")) as FqLintResult

    topic:
    tuple(task.process, 'fq', eval("fq lint --version | sed 's/fq-lint //; s/ .*//'")) >> 'versions'

    script:
    def args   = task.ext.args ?: sample.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    fq lint \\
        $args \\
        ${sample.reads.join(' ')} > ${prefix}.fq_lint.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.fq_lint.txt
    """
}
