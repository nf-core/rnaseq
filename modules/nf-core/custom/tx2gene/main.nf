process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    meta        : Map
    gtf         : Path
    quants      : List<Path>
    quant_type  : String
    id          : String
    extra       : String

    stage:
    stageAs "quants/*", quants

    output:
    file("*tx2gene.tsv")

    script:
    template 'tx2gene.py'

    stub:
    """
    touch ${meta.id}.tx2gene.tsv
    """
}
