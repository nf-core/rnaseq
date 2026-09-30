nextflow.enable.types = true

record MultiqcWriteFileInput {
    id:      String
    meta:    Map
    name:    String
    content: String
}

process MULTIQC_WRITE_FILE {
    tag "${sample.id}"

    input:
    sample: MultiqcWriteFileInput

    output:
    record(id: sample.id, meta: sample.meta, file: file(sample.name))

    exec:
    task.workDir.resolve(sample.name) << sample.content
}
