// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process MULTIQC_CUSTOM_FAIL_MAPPED {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    val  fail_mapped
    path header
    val  options

    output:
    path "*.tsv"

    script:
    if (fail_mapped.size() > 0) {
        """
        sed 's/MIN_MAPPED_READS/${params.min_mapped_reads}/g' $header > header.tmp.txt
        echo "${fail_mapped.join('\n')}" > fail.tmp.txt
        cat header.tmp.txt fail.tmp.txt  > fail_mapping_mqc.tsv
        """
    } else {
        """
        touch fail_mapping_mqc.tsv
        """
    }
}
