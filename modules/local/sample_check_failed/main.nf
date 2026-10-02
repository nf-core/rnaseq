process SAMPLE_CHECK_FAILED {
    tag "${id}"

    input:
    tuple val(id), val(message)

    exec:
    throw new nextflow.exception.ProcessFailedException(message)
}
