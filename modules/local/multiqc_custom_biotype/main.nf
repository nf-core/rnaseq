process MULTIQC_CUSTOM_BIOTYPE {
    tag "$meta.id"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    meta    : Map
    count   : Path
    header  : Path

    output:
    file("*.tsv")

    topic:
    file('versions.yml') >> 'versions'
    file("*.tsv") >> 'logs'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,7 $count | tail -n +3 | cat $header - >> ${prefix}.biotype_counts_mqc.tsv

    mqc_features_stat.py \\
        ${prefix}.biotype_counts_mqc.tsv \\
        -s $meta.id \\
        -f rRNA \\
        -o ${prefix}.biotype_counts_rrna_mqc.tsv
    """
}
