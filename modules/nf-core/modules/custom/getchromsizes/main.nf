process CUSTOM_GETCHROMSIZES {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15--h1170115_1' :
        'quay.io/biocontainers/samtools:1.15--h1170115_1' }"

    input:
    path fasta

    output:
    path '*.sizes'      , emit: sizes
    path '*.fai'        , emit: fai
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        custom: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
