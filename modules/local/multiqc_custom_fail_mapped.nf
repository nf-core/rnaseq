// Import generic module functions
include { saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_FAIL_MAPPED {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val fail_mapped

    output:
    path "*.tsv"

    script:
    if (fail_mapped.size() > 0) {
        """
        echo "Sample\tSTAR uniquely mapped reads (%)" > fail_mapped_samples_mqc.tsv
        echo "${fail_mapped.join('\n')}" >> fail_mapped_samples_mqc.tsv
        """
    } else {
        """
        touch fail_mapped_samples_mqc.tsv
        """
    }
}
