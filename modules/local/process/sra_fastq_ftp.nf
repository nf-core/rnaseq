// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Download SRA data via FTP
 */
process SRA_FASTQ_FTP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }
    
    input:
    tuple val(meta), val(fastq)

    output:
    tuple val(meta), path("*fastq.gz"), emit: fastq
    tuple val(meta), path("*md5")     , emit: md5

    script:    
    if (meta.single_end) {
        """
        curl -L ${fastq[0]} -o ${meta.id}.fastq.gz
        echo "${meta.md5_1} ${meta.id}.fastq.gz" > ${meta.id}.fastq.gz.md5
        md5sum -c ${meta.id}.fastq.gz.md5
        """
    } else {
        """
        curl -L ${fastq[0]} -o ${meta.id}_1.fastq.gz
        echo "${meta.md5_1} ${meta.id}_1.fastq.gz" > ${meta.id}_1.fastq.gz.md5
        md5sum -c ${meta.id}_1.fastq.gz.md5

        curl -L ${fastq[1]} -o ${meta.id}_2.fastq.gz
        echo "${meta.md5_2} ${meta.id}_2.fastq.gz" > ${meta.id}_2.fastq.gz.md5
        md5sum -c ${meta.id}_2.fastq.gz.md5
        """
    }
}
