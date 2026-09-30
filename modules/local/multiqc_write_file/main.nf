nextflow.enable.types = true

include { MultiqcWriteFileInput } from '../../nf-core/types'

process MULTIQC_WRITE_FILE {
    tag "${sample.id}"

    input:
    sample: MultiqcWriteFileInput

    output:
    record(id: sample.id, meta: sample.meta, file: file(sample.name))

    exec:
    task.workDir.resolve(sample.name) << sample.content
}
