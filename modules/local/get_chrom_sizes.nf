// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process GET_CHROM_SIZES {
    tag "$fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    path fasta

    output:
    path '*.sizes'      , emit: sizes
    path '*.fai'        , emit: fai
    path  "versions.yml", emit: version

    script:
    def software = 'samtools'
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
