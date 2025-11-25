process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data':
        'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01' }"

    input:
    tuple val(meta), path(kraken_report)
    path database

    output:
    tuple val(meta), path(bracken_report)        , emit: reports
    tuple val(meta), path(bracken_kraken_style_report), emit: txt
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${kraken_report}' \\
        -o '${bracken_report}' \\
        -w '${bracken_kraken_style_report}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.tsv"
    bracken_kraken_style_report = "${prefix}.kraken2.report_bracken.txt"
    """
    touch ${prefix}.tsv
    touch ${bracken_kraken_style_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """
}
