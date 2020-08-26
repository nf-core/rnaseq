// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process SALMON_TX2GENE {
    tag "$gtf"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path ("salmon/*")
    path gtf
    val options

    output:
    path "*.csv", emit: csv

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tx2gene.py \\
        --gtf $gtf \\
        --salmon salmon \\
        --id $params.fc_group_features \\
        --extra $params.fc_extra_attributes \\
        -o salmon_tx2gene.csv
    """
}
