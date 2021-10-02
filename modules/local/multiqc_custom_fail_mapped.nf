// Import generic module functions
include { saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_FAIL_MAPPED {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    memory 100.MB

    input:
    val fail_mapped

    output:
    path "*.tsv", emit: tsv

    exec:
    def header = [
        "Sample",
        "STAR uniquely mapped reads (%)"
    ]
    def contents = "${header.join('\t')}\n"
    def mqc_file = task.workDir.resolve("fail_mapped_samples_mqc.tsv")
    if (fail_mapped.size() > 0) {
        contents += fail_mapped.join('\n')
        mqc_file.text = contents
    }
}
