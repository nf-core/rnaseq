nextflow.enable.types = true

include { FastaInput } from '../../types'

record KallistoIndexResult {
    id:    String
    meta:  Map
    index: Path
}

process KALLISTO_INDEX {
    tag "$sample.fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e21d2cff2526b0995996977c057f0c17844073781f86a18b15fef178a97ed7cb/data':
        'community.wave.seqera.io/library/kallisto:0.52.0--31c771060d82d25c' }"

    input:
    sample: FastaInput

    output:
    record(id: sample.id, meta: sample.meta, index: file('kallisto')) as KallistoIndexResult

    topic:
    tuple(task.process, 'kallisto', eval('kallisto 2>&1 | head -1 | sed "s/^kallisto //; s/Usage.*//"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    kallisto \\
        index \\
        $args \\
        -i kallisto \\
        $sample.fasta
    """

    stub:
    """
    mkdir kallisto
    """
}
