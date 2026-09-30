nextflow.enable.types = true

record MultiqcConcatenateTablesInput {
    id:    String
    meta:  Map
    name:  String
    skip:  Integer
    files: List<Path>
}

record MultiqcConcatenateTablesResult {
    id:   String
    meta: Map
    file: Path
}

process MULTIQC_CONCATENATE_TABLES {
    tag "${sample.id}"
    label 'process_single'

    input:
    sample: MultiqcConcatenateTablesInput

    output:
    record(id: sample.id, meta: sample.meta, file: file(sample.name)) as MultiqcConcatenateTablesResult

    script:
    def first = sample.files.first().name
    def rest  = sample.files.tail().collect { f -> f.name }.join(' ')
    """
    cat ${first} > ${sample.name}
    for f in ${rest}; do
        tail -n +${sample.skip + 1} "\$f" >> ${sample.name}
    done
    """

    stub:
    """
    touch ${sample.name}
    """
}
