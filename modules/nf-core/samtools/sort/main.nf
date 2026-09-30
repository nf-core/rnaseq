nextflow.enable.types = true

process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple(meta: Map, bam: List<Path>)
    tuple(meta2: Map, fasta: Path?, fai: Path?)
    index_format: String

    stage:
    stageAs bam, '?/*'

    output:
    record(
        meta: meta,
        bam:  file("${task.ext.prefix ?: meta.id}.bam", optional: true),
        cram: file("${task.ext.prefix ?: meta.id}.cram", optional: true),
        sam:  file("${task.ext.prefix ?: meta.id}.sam", optional: true),
        bai:  file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.bai", optional: true),
        csi:  file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.csi", optional: true),
        crai: file("${task.ext.prefix ?: meta.id}.{bam,cram,sam}.crai", optional: true)
    )

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam")
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
    // A lone file arrives as a Path, which iterates over its name components.
    def raw_bam = bam as Object
    def bams = raw_bam instanceof List ? (raw_bam as List<Path>) : [raw_bam as Path]
    def is_sam = bams[0].name.endsWith('.sam')
    if (index_format) {
        if (!(index_format in ['bai', 'csi', 'crai'])) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }
    if (bams.join(' ') == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    def input_source = is_sam ? bams.join(' ') : "-"
    def pre_command = is_sam ? "" : "samtools cat ${bams.join(' ')} | "

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam")
        ? "sam"
        : args.contains("--output-fmt cram")
            ? "cram"
            : "bam"

    if (index_format) {
        if (!(index_format in ['bai', 'csi', 'crai'])) {
            error("Index format not one of bai, csi, crai.")
        }
        else if (extension == "sam") {
            error("Indexing not compatible with SAM output")
        }
    }

    def index = index_format ? "touch ${prefix}.${extension}.${index_format}" : ""

    """
    touch ${prefix}.${extension}
    ${index}
    """
}
