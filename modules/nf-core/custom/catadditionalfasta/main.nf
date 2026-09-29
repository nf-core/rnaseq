nextflow.enable.types = true

process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple(meta: Map, fasta: Path, gtf: Path)
    tuple(meta2: Map, add_fasta: Path)
    biotype: String

    output:
    record(
        meta:  meta,
        fasta: file("out/${task.ext.prefix ?: meta.id}.fasta"),
        gtf:   file("out/${task.ext.prefix ?: meta.id}.gtf")
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'fasta2gtf.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/${prefix}.fasta
    touch out/${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
