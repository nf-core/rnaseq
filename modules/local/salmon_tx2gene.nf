// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

process SALMON_TX2GENE {
    tag "$gtf"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path ("salmon/*")
    path gtf

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tx2gene.py \\
        --gtf $gtf \\
        --salmon salmon \\
        --id $params.gtf_group_features \\
        --extra $params.gtf_extra_attributes \\
        -o salmon_tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
