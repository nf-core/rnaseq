nextflow.enable.types = true

process SAMTOOLS_VIEW {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple(meta: Map, input: Path, index: Path?)
    tuple(meta2: Map, fasta: Path?, fai: Path?)
    tuple(meta3: Map, qname: Path?)
    tuple(meta4: Map, bed: Path?)
    index_format: String

    output:
    record(
        meta:             meta,
        bam:              file("${task.ext.prefix ?: meta.id}.bam", optional: true),
        cram:             file("${task.ext.prefix ?: meta.id}.cram", optional: true),
        sam:              file("${task.ext.prefix ?: meta.id}.sam", optional: true),
        bai:              file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.bai", optional: true),
        csi:              file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.csi", optional: true),
        crai:             file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.crai", optional: true),
        unselected:       file("${task.ext.prefix ?: meta.id}.unselected.{bam,cram,sam}", optional: true),
        unselected_index: file("${task.ext.prefix ?: meta.id}.unselected.{bam,cram,sam}.{csi,crai}", optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def file_type = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : input.getExtension()

    def output_file = index_format ? "${prefix}.${file_type}##idx##${prefix}.${file_type}.${index_format} --write-index" : "${prefix}.${file_type}"
    // Can't choose index type of unselected file
    def readnames = qname ? "--qname-file ${qname} --output-unselected ${prefix}.unselected.${file_type}" : ""
    def bedfile = bed ? "-L ${bed}" : ""

    if ("${input}" == "${prefix}.${file_type}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    if (index_format) {
        if (!(index_format in ['bai', 'csi', 'crai'])) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (file_type == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        view \\
        --threads ${task.cpus - 1} \\
        ${reference} \\
        ${readnames} \\
        ${bedfile} \\
        ${args} \\
        -o ${output_file} \\
        ${input} \\
        ${args2}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt bam")
            ? "bam"
            : args.contains("--output-fmt cram")
                ? "cram"
                : input.getExtension()
    def default_index_format = file_type == "bam"
        ? "csi"
        : file_type == "cram" ? "crai" : ""
    def index_command = index_format ? "touch ${prefix}.${file_type}.${index_format}" : args.contains("--write-index") ? "touch ${prefix}.${file_type}.${default_index_format}" : ""
    def unselected = qname ? "touch ${prefix}.unselected.${file_type}" : ""
    // Can't choose index type of unselected file
    def unselected_index = qname && (args.contains("--write-index") || index_format) ? "touch ${prefix}.unselected.${file_type}.${default_index_format}" : ""

    if ("${input}" == "${prefix}.${file_type}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    if (index_format) {
        if (!(index_format in ['bai', 'csi', 'crai'])) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (file_type == "sam") {
            error("Indexing not compatible with SAM output.")
        }
    }
    """
    touch ${prefix}.${file_type}
    ${index_command}
    ${unselected}
    ${unselected_index}
    """
}
