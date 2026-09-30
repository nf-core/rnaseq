nextflow.enable.types = true

include { SeqkitReplaceInput } from '../../types'

process SEQKIT_REPLACE {
    tag "${sample.meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    sample: SeqkitReplaceInput
    out_ext: String

    output:
    record(id: sample.id, meta: sample.meta, fastx: file("*.fast*"))

    topic:
    tuple(task.process, 'seqkit', eval("seqkit version | sed 's/^.*v//'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def extension = "fastq"
    if ("${sample.fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.faa|.+\.faa.gz/) {
        extension = "fasta"
    }
    def isgz = ""
    if ("${sample.fastx}" ==~ /.+\.gz/) {
        isgz = ".gz"
    }
    def endswith = out_ext ?: "${extension}${isgz}"
    """
    seqkit \\
        replace \\
        ${args} \\
        --threads ${task.cpus} \\
        -i ${sample.fastx} \\
        -o ${prefix}.${endswith}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def extension = "fastq"
    if ("${sample.fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.faa|.+\.faa.gz/) {
        extension = "fasta"
    }
    def isgz = ""
    if ("${sample.fastx}" ==~ /.+\.gz/) {
        isgz = ".gz"
    }
    def endswith = out_ext ?: "${extension}${isgz}"

    def create_cmd = endswith.endsWith('gz') ? "echo '' | gzip >" : "touch"
    """
    echo ${args}

    ${create_cmd} ${prefix}.${endswith}
    """
}
