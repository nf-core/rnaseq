// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52d3fbb86c48b8c39a989733ca97bfd84cf665a9ab279be7cfbbb36f5a53c46a/data'
        : 'community.wave.seqera.io/library/bedtools:2.31.1--897920cd303d6413' }"

    input:
    tuple val(meta), path(intervals), val(scale)
    path sizes
    val extension
    val sort

    output:
    tuple val(meta), path("*.${extension}"), emit: genomecov
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

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
