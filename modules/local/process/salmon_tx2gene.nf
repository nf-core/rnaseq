// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process SALMON_TX2GENE {
    tag "$gtf"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity') {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path ("salmon/*")
    path gtf
    
    output:
    path "*.tsv"

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tx2gene.py \\
        --gtf $gtf \\
        --salmon salmon \\
        --id $params.fc_group_features \\
        --extra $params.fc_extra_attributes \\
        -o salmon_tx2gene.tsv
    """
}
