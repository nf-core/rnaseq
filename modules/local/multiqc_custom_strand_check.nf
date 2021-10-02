// Import generic module functions
include { saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_STRAND_CHECK {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    memory 100.MB

    input:
    val fail_strand

    output:
    path "*.tsv", emit: tsv

    exec:
    // Generate file contents
    def header = [
        "Sample",
        "Provided strandedness",
        "Inferred strandedness",
        "Sense (%)",
        "Antisense (%)",
        "Undetermined (%)"
    ]
    def contents = ''
    if (fail_strand.size() > 0) {
        contents += "${header.join('\t')}\n"
        contents += fail_strand.join('\n')
    }

    // Write to file
    def mqc_file = task.workDir.resolve("fail_strand_check_mqc.tsv")
    mqc_file.text = contents
}
