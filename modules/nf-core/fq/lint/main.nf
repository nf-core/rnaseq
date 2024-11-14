process FQ_LINT {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'terminate'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fq_lint.txt"), emit: lint
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq lint \\
        $args \\
        $fastq > ${prefix}.fq_lint.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq lint --version | sed 's/fq-lint //g'))
    END_VERSIONS

    if ! tail -n 1 ${prefix}.fq_lint.txt | grep -q 'fq-lint end'; then
        echo "ERROR: Linting failure detected for ${meta.id}. See ${prefix}.fq_lint.txt for details."
        exit 1
    fi
    """
}
