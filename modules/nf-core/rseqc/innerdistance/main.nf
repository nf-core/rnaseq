process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    meta    : Map
    bam     : Path
    bed     : Path

    output:
    distance    : Path? = file("*distance.txt")
    freq        : Path? = file("*freq.txt")
    mean        : Path? = file("*mean.txt")
    pdf         : Path? = file("*.pdf")
    rscript     : Path? = file("*.r")

    topic:
    file('versions.yml') >> 'versions'
    file('*freq.txt') >> 'logs'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!meta.single_end) {
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            $args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt
        """
    } else {
        """
        """
    }
}
