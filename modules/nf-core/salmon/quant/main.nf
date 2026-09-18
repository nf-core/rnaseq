process SALMON_QUANT {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1ce42a19f9e7135babf14432e80b33ec717a13f14734d150d53347e979919629/data'
        : 'community.wave.seqera.io/library/salmon:2.7.0--74784226202c61b9'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index), path(gtf), path(transcript_fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional: true
    tuple val(meta), path("*lib_format_counts.json"), emit: lib_format_counts, optional: true
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | sed -e "s/salmon //g"'), topic: versions, emit: versions_salmon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // salmon's -a takes a BAM of reads already aligned to the transcriptome; anything else is reads mode
    def alignment_mode = "${reads instanceof List ? reads[0] : reads}".endsWith('.bam')

    def transcript_fasta_file = transcript_fasta instanceof List ? (transcript_fasta ? transcript_fasta[0] : null) : transcript_fasta
    def index_dir = index instanceof List ? (index ? index[0] : null) : index

    def reference
    def input_reads
    if (alignment_mode) {
        if (!transcript_fasta_file) {
            error("[Salmon Quant] Alignment mode needs 'transcript_fasta' to be provided as the reference (BAM input detected for sample '${meta.id}').")
        }
        reference = "-t ${transcript_fasta}"
        input_reads = "-a ${reads}"
    }
    else {
        if (!index_dir) {
            error("[Salmon Quant] Reads mode needs 'index' to be provided as a salmon index directory (no BAM input detected for sample '${meta.id}').")
        }
        def reads1 = []
        def reads2 = []
        meta.single_end ? [reads].flatten().each { r -> reads1 << r } : reads.eachWithIndex { v, ix -> (ix & 1 ? reads2 : reads1) << v }
        reference = "--index ${index}"
        input_reads = meta.single_end ? "-r ${reads1.join(" ")}" : "-1 ${reads1.join(" ")} -2 ${reads2.join(" ")}"
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
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}_meta_info.json
    touch ${prefix}_lib_format_counts.json
    """
}
