nextflow.enable.types = true

process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6397750e9730a3fbcc5b4c43f14bd141c64c723fd7dad80e47921a68a7c3cd21/data'
        : 'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b'}"

    input:
    record(id: String, meta: Map, intervals: Path, scale: Float)
    sizes: Path?
    extension: String
    sort: Boolean

    output:
    record(id: id, meta: meta, bedgraph: file("*.${extension}"))

    topic:
    tuple(task.process, 'bedtools', eval("bedtools --version | sed -e 's/bedtools v//g'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()
    args += scale > 0 && scale != 1 ? " -scale ${scale}" : ""
    if (!args_list.contains('-bg') && (scale > 0 && scale != 1)) {
        args += " -bg"
    }
    // Sorts output file by chromosome and position using additional options for performance and consistency
    // See https://www.biostars.org/p/66927/ for further details
    def buffer = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def sort_cmd = sort ? "| LC_ALL=C sort --parallel=${task.cpus} ${buffer} -k1,1 -k2,2n" : ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    if (intervals.name =~ /\.bam/) {
        """
        bedtools \\
            genomecov \\
            -ibam ${intervals} \\
            ${args} \\
            ${sort_cmd} \\
            > ${prefix}.${extension}
        """
    }
    else {
        """
        bedtools \\
            genomecov \\
            -i ${intervals} \\
            -g ${sizes} \\
            ${args} \\
            ${sort_cmd} \\
            > ${prefix}.${extension}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.${extension}
    """
}
