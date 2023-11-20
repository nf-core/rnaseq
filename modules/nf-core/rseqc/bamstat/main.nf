process RSEQC_BAMSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt
    tuple val("${task.process}"), val('rseqc'), cmd("bam_stat.py --version | sed -e 's/bam_stat.py //g'"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_stat.py \\
        -i $bam \\
        $args \\
        > ${prefix}.bam_stat.txt
    """
}
