process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta), path(fai)
    val index_format

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam, optional: true
    tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
    tuple val(meta), path("${prefix}.sam"), emit: sam, optional: true
    tuple val(meta), path("${prefix}.${extension}.{crai,csi,bai}"), emit: index, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt cram")
            ? "cram"
            : "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    //setting default values
    def write_index = ""
    def output_file = "${prefix}.${extension}"

    // Update if index is requested
    if (index_format != '' && index_format) {
        write_index = "--write-index"
        output_file = "${prefix}.${extension}##idx##${prefix}.${extension}.${index_format}"
    }
    def is_sam = (bam instanceof List ? bam[0] : bam).name.endsWith('.sam')
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    def input_source = is_sam ? "${bam}" : "-"
    def pre_command = is_sam ? "" : "samtools cat ${bam} | "

    """
    ${pre_command}samtools sort \\
        ${args} \\
        -T ${prefix} \\
        --threads ${task.cpus} \\
        ${reference} \\
        -o ${output_file} \\
        ${write_index} \\
        ${input_source}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt cram")
            ? "cram"
            : "bam"

    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }

    index = index_format ? "touch ${prefix}.${extension}.${index_format}" : ""

    """
    touch ${prefix}.${extension}
    ${index}
    """
}
