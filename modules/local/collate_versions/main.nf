nextflow.enable.types = true

process COLLATE_VERSIONS {

    input:
    name: String
    entries: List<String>

    output:
    file(name)

    exec:
    def collated = task.workDir.resolve(name)
    entries.each { entry ->
        collated << entry
        collated << '\n'
    }
}
