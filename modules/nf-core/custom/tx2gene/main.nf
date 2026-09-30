nextflow.enable.types = true

process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    record(id: String, meta: Map, quants: List<Path>)
    gtf: Path
    quant_type: String
    gtf_id_attribute: String
    extra: String?

    stage:
    stageAs quants, 'quants/*'

    output:
    record(id: id, meta: meta, tx2gene: file("*tx2gene.tsv"))

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'tx2gene.py'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
