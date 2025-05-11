process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    meta    : Map
    bam     : Path

    output:
    seq_xls : Path = file("*seq.DupRate.xls")
    pos_xls : Path = file("*pos.DupRate.xls")
    pdf     : Path = file("*.pdf")
    rscript : Path = file("*.r")

    topic:
    file('versions.yml') >> 'versions'
    file('*pos.DupRate.xls') >> 'logs'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
        $args
    """
}
