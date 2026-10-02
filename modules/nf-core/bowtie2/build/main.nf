nextflow.enable.types = true

include { FastaInput } from '../../types'

record Bowtie2BuildResult {
    id:    String
    meta:  Map
    index: Path
}

process BOWTIE2_BUILD {
    tag "${sample.fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    sample: FastaInput

    output:
    record(id: sample.id, meta: sample.meta, index: file('bowtie2')) as Bowtie2BuildResult

    topic:
    tuple(task.process, 'bowtie2', eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build $args --threads $task.cpus ${sample.fasta} bowtie2/${sample.fasta.baseName}
    """

    stub:
    """
    mkdir bowtie2
    touch bowtie2/${sample.fasta.baseName}.{1..4}.bt2
    touch bowtie2/${sample.fasta.baseName}.rev.{1,2}.bt2
    """
}
