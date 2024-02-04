process RSEQC_TIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::rseqc=5.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.xls"), emit: xls
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tin.py \\
        -i $bam \\
        -r $bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(tin.py --version | sed -e "s/tin.py //g")
    END_VERSIONS
    """
}
