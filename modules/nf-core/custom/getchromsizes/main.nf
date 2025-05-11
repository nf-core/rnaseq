process CUSTOM_GETCHROMSIZES {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    fasta   : Path

    output:
    sizes   : Path = file("*.sizes")
    fai     : Path = file("*.fai")
    gzi     : Path? = file("*.gzi")

    script:
    def args = task.ext.args ?: ''
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """

    stub:
    """
    touch ${fasta}.fai
    touch ${fasta}.sizes
    """
}
