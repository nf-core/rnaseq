// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/umi_tools:1.1.1--py38h0213d0e_1"
    } else {
        container "quay.io/biocontainers/umi_tools:1.1.1--py38h0213d0e_1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $options.args \\
            > ${prefix}.umi_extract.log

        umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
        """
    }  else {
        """
        umi_tools \\
            extract \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.umi_extract_1.fastq.gz \\
            --read2-out=${prefix}.umi_extract_2.fastq.gz \\
            $options.args \\
            > ${prefix}.umi_extract.log

        umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
        """
    }
}
