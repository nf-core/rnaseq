// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process MULTIQC_CUSTOM_STRAND_CHECK {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    val fail_strand
    val options

    output:
    path "*.tsv"

    script:
    if (fail_strand.size() > 0) {
        """
        echo "Sample\tProvided strandedness\tInferred strandedness\tSense (%)\tAntisense (%)\tUndetermined (%)" > fail_strand_check_mqc.tsv
        echo "${fail_strand.join('\n')}" >> fail_strand_check_mqc.tsv
        """
    } else {
        """
        touch fail_strand_check_mqc.tsv
        """
    }
}
