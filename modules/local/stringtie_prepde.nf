// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRINGTIE_PREPDE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"
    }

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gene_count_matrix.csv")      , emit: counts_gene
    tuple val(meta), path("*.transcript_count_matrix.csv"), emit: counts_transcript
    path  "*.version.txt"                                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "${meta.id} ${gtf}" > sample_lst.txt

    prepDE.py \\
        -i sample_lst.txt \\
        -g ${prefix}.gene_count_matrix.csv \\
        -t ${prefix}.transcript_count_matrix.csv \\
        -l $meta.read_length \\
        $options.args

    echo \$(stringtie --version 2>&1) > ${software}.version.txt
    """
}
