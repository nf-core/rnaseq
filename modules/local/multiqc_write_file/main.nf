nextflow.enable.types = true

process MULTIQC_WRITE_FILE {
    tag "${id}"

    input:
    record(id: String, meta: Map, name: String, content: String)

    output:
    record(id: id, meta: meta, file: file(name))

    exec:
    task.workDir.resolve(name) << content
}
