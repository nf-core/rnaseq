nextflow.enable.types = true

record SamtoolsFaidxInput {
    id:    String
    meta:  Map
    fasta: Path
    fai:   Path?
}

process SAMTOOLS_FAIDX {
    tag "${sample.fasta}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    sample: SamtoolsFaidxInput
    get_sizes: Boolean

    output:
    record(
        id:    sample.id,
        meta:  sample.meta,
        fa:    file('*.{fa,fasta}', optional: true),
        sizes: file('*.sizes', optional: true),
        fai:   file('*.fai', optional: true),
        gzi:   file('*.gzi', optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def get_sizes_command = get_sizes ? "cut -f 1,2 ${sample.fasta}.fai > ${sample.fasta}.sizes" : ''
    """
    samtools \\
        faidx \\
        ${sample.fasta} \\
        ${args}

    ${get_sizes_command}
    """

    stub:
    def match = (task.ext.args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    def get_sizes_command = get_sizes ? "touch ${sample.fasta}.sizes" : ''
    """
    ${fastacmd}
    touch ${sample.fasta}.fai
    if [[ "${sample.fasta.extension}" == "gz" ]]; then
        touch ${sample.fasta}.gzi
    fi

    ${get_sizes_command}
    """
}
