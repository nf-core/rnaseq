process GTF_FILTER {
    tag "$fasta"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    fasta   : Path
    gtf     : Path

    output:
    filtered_gtf: Path = file("*.filtered.gtf")

    script: // filter_gtf.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        --prefix ${fasta.baseName}
    """
}
