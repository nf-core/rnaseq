process CUSTOM_MULTIQCCUSTOMBIOTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(count)
    tuple val(meta2), path(header)

    output:
    tuple val(meta), path("*biotype_counts_mqc.tsv")      , emit: tsv
    tuple val(meta), path("*biotype_counts_rrna_mqc.tsv") , emit: rrna
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'mqc_features_stat.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.biotype_counts_mqc.tsv
    touch ${prefix}.biotype_counts_rrna_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
