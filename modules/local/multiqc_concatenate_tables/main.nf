nextflow.enable.types = true

process MULTIQC_CONCATENATE_TABLES {
    tag "${id}"

    input:
    record(id: String, meta: Map, name: String, skip: Integer, files: List<Path>)

    output:
    record(id: id, meta: meta, file: file(name))

    exec:
    def combined = task.workDir.resolve(name)
    combined << files.first().text
    files.tail().each { f ->
        def lines = f.readLines()
        combined << lines.subList(skip, lines.size()).collect { line -> "${line}\n" }.join('')
    }
}
