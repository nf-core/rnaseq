// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/umi_tools:1.0.1--py37h516909a_1"
    container "https://depot.galaxyproject.org/singularity/umi_tools:1.0.1--py37h516909a_1"

    conda (params.conda ? "bioconda::umi_tools=1.0.1" : null)

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.rm_umi.fastq.gz \\
            $ioptions.args \\
            > ${prefix}.umi.log

        umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
        """
    }  else {
        """
        umi_tools \\
            extract \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.rm_umi_1.fastq.gz \\
            --read2-out=${prefix}.rm_umi_2.fastq.gz \\
            $ioptions.args \\
            > ${prefix}.umi.log

        umi_tools --version | sed -e "s/UMI-tools version: //g" > ${software}.version.txt
        """
    }
}
