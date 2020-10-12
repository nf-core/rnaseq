// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process MULTIQC_CUSTOM_BIOTYPE {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/python:3.8.3"
    //container  https://depot.galaxyproject.org/singularity/python:3.8.3

    conda (params.conda ? "conda-forge::python=3.8.3" : null)

    input:
    tuple val(meta), path(count)
    path  header
    val   options

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    cut -f 1,7 $count | tail -n +3 | cat $header - >> ${prefix}.biotype_counts_mqc.tsv
    mqc_features_stat.py ${prefix}.biotype_counts_mqc.tsv -s $meta.id -f rRNA -o ${prefix}.biotype_counts_rrna_mqc.tsv
    """
}
