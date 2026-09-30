nextflow.enable.types = true

record MultiqcConcatenateTablesInput {
    id:    String
    meta:  Map
    name:  String
    skip:  Integer
    files: List<Path>
}

process MULTIQC_CONCATENATE_TABLES {
    tag "${sample.id}"

    input:
    sample: MultiqcConcatenateTablesInput

    output:
    record(id: sample.id, meta: sample.meta, file: file(sample.name))

    exec:
    def combined = task.workDir.resolve(sample.name)
    combined << sample.files.first().text
    sample.files.tail().each { f ->
        def lines = f.readLines()
        combined << lines.subList(sample.skip, lines.size()).collect { line -> "${line}\n" }.join('')
    }
}
