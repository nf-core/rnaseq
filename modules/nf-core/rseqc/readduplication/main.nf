process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*seq.DupRate.xls"), emit: seq_xls
    tuple val(meta), path("*pos.DupRate.xls"), emit: pos_xls
    tuple val(meta), path("*.pdf")           , emit: pdf
    tuple val(meta), path("*.r")             , emit: rscript
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_duplication.py --version | sed -e "s/read_duplication.py //g")
    END_VERSIONS
    """
}
