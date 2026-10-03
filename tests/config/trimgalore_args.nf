process TRIMGALORE {
    input:
    val trigger

    output:
    val args, emit: args

    script:
    args = task.ext.args
    """
    true
    """
}

workflow FASTQ_FASTQC_UMITOOLS_TRIMGALORE {
    take:
    trigger

    main:
    TRIMGALORE(trigger)

    emit:
    TRIMGALORE.out.args
}

workflow TRIMGALORE_CONFIG_TEST {
    take:
    trigger

    main:
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE(trigger)

    emit:
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out[0]
}
