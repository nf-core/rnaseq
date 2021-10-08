// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

process GTF_GENE_FILTER {
    tag "$fasta"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path fasta
    path gtf

    output:
    path "*.gtf"       , emit: gtf
    path "versions.yml", emit: versions

    script: // filter_gtf_for_genes_in_genome.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf_for_genes_in_genome.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        -o ${fasta.baseName}_genes.gtf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
