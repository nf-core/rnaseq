process TX2GENE {
    tag "$gtf"
    label "process_low"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path ("quants/*")
    val quant_type
    path gtf

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    tx2gene.py \\
        --quant_type $quant_type \\
        --gtf $gtf \\
        --quants quants \\
        --id $params.gtf_group_features \\
        --extra $params.gtf_extra_attributes \\
        -o tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
