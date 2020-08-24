// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/fastqc:0.11.9--0"
    //container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"

    conda (params.conda ? "bioconda::fastqc=0.11.9" : null)

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "*.version.txt", emit: version

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}.${ioptions.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        fastqc $ioptions.args --threads $task.cpus ${prefix}.fastq.gz
        fastqc --version | sed -e "s/FastQC v//g" > ${software}.version.txt
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        fastqc $ioptions.args --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
        fastqc --version | sed -e "s/FastQC v//g" > ${software}.version.txt
        """
    }
}
