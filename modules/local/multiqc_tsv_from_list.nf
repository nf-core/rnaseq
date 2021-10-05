// Import generic module functions
include { saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

process MULTIQC_TSV_FROM_LIST {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    memory 100.MB

    input:
    val tsv_data   // [ ['foo', 1], ['bar', 1] ]
    val header     // [ 'name', 'number' ]
    val out_prefix

    output:
    path "*.tsv"

    exec:
    // Generate file contents
    def contents = ""
    if (tsv_data.size() > 0) {
        contents += "${header.join('\t')}\n"
        contents += tsv_data.join('\n')
    }

    // Write to file
    def mqc_file = task.workDir.resolve("${out_prefix}_mqc.tsv")
    mqc_file.text = contents
}
