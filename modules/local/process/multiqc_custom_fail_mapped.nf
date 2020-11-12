// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_FAIL_MAPPED {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "biocontainers/biocontainers:v1.2.0_cv1"
    
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
