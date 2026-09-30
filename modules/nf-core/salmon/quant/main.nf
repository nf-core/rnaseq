nextflow.enable.types = true

include { ReadsInput } from '../../types'

record SalmonQuantSample {
    id:                String
    meta:              Map
    quant_dir:         Path
    json_info:         Path?
    lib_format_counts: Path?
}

process SALMON_QUANT {
    tag "${sample.meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1ce42a19f9e7135babf14432e80b33ec717a13f14734d150d53347e979919629/data'
        : 'community.wave.seqera.io/library/salmon:2.7.0--74784226202c61b9'}"

    input:
    sample: ReadsInput
    index: Path?
    gtf: Path
    transcript_fasta: Path?

    output:
    record(
        id:                sample.id,
        meta:              sample.meta,
        quant_dir:         file("${prefix}"),
        json_info:         file("*info.json", optional: true),
        lib_format_counts: file("*lib_format_counts.json", optional: true)
    ) as SalmonQuantSample

    topic:
    tuple(task.process, 'salmon', eval('salmon --version | sed -e "s/salmon //g"')) >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${sample.meta.id}"

    // salmon's -a takes a BAM of reads already aligned to the transcriptome; anything else is reads mode
    // A lone file arrives as a Path, which iterates over its name components.
    def raw_reads = sample.reads as Object
    def reads_list = raw_reads instanceof List ? (raw_reads as List<Path>) : [raw_reads as Path]
    def alignment_mode = "${reads_list[0]}".endsWith('.bam')

    def transcript_fasta_file = transcript_fasta
    def index_dir = index

    def reference
    def input_reads
    if (alignment_mode) {
        if (!transcript_fasta_file) {
            error("[Salmon Quant] Alignment mode needs 'transcript_fasta' to be provided as the reference (BAM input detected for sample '${sample.meta.id}').")
        }
        reference = "-t ${transcript_fasta}"
        input_reads = "-a ${reads_list.join(' ')}"
    }
    else {
        if (!index_dir) {
            error("[Salmon Quant] Reads mode needs 'index' to be provided as a salmon index directory (no BAM input detected for sample '${sample.meta.id}').")
        }
        def reads1 = sample.meta.single_end ? reads_list.join(' ') : reads_list.findAll { r -> reads_list.indexOf(r) % 2 == 0 }.join(' ')
        def reads2 = reads_list.findAll { r -> reads_list.indexOf(r) % 2 == 1 }.join(' ')
        reference = "--index ${index}"
        input_reads = sample.meta.single_end ? "-r ${reads1}" : "-1 ${reads1} -2 ${reads2}"
    }

    """
    salmon quant \\
        --geneMap ${gtf} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input_reads} \\
        ${args} \\
        -o ${prefix}

    if [ -f ${prefix}/aux_info/meta_info.json ]; then
        cp ${prefix}/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi
    if [ -f ${prefix}/lib_format_counts.json ]; then
        cp ${prefix}/lib_format_counts.json "${prefix}_lib_format_counts.json"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}_meta_info.json
    touch ${prefix}_lib_format_counts.json
    """
}
