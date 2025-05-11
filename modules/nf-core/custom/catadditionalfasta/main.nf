process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    fasta       : Path
    gtf         : Path
    add_fasta   : Path
    biotype     : String
    prefix      : String = ''

    output:
    fasta       : Path = file("*/*.fasta")
    gtf         : Path = file("*/*.gtf")

    script:
    template 'fasta2gtf.py'
}
