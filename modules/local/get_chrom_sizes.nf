// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

process GET_CHROM_SIZES {
    tag "$fasta"

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    path fasta

    output:
    path '*.sizes'     , emit: sizes
    path '*.fai'       , emit: fai
    path "versions.yml", emit: versions

    script:
    """
    samtools \\
        faidx \\
        $fasta

    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
