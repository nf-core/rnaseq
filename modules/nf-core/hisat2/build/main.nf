process HISAT2_BUILD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a16b8102041e76f25358477471bafe813828313f7a8787066f390354cd7b8b7c/data'
        : 'community.wave.seqera.io/library/hisat2:2.2.3--2616fa83d3b9d8f8'}"

    input:
    tuple val(meta), path(fasta), path(gtf), path(splicesites)
    val hisat2_memory_input

    output:
    tuple val(meta), path("hisat2"), emit: index
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n 's/.*version \\([^ ]*\\).*/\\1/p'"), emit: versions_hisat2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (!task.memory) {
        error("[HISAT2 index build] No memory specified for process. Please configure memory for 'process_high' label.")
    }
    def avail_mem = task.memory.toGiga()
    def hisat2_build_memory = hisat2_memory_input ? (hisat2_memory_input as MemoryUnit).toGiga() : Integer.MAX_VALUE

    def ss = ''
    def exon = ''
    def extract_exons = ''

    if (avail_mem >= hisat2_build_memory) {
        log.info("[HISAT2 index build] ${avail_mem} GB available, using splice sites and exons to build HISAT2 index")
        extract_exons = gtf ? "hisat2_extract_exons.py ${gtf} > ${gtf.baseName}.exons.txt" : ""
        ss = splicesites ? "--ss ${splicesites}" : ""
        exon = gtf ? "--exon ${gtf.baseName}.exons.txt" : ""
    }
    else {
        log.info("[HISAT2 index build] Only ${avail_mem} GB available (< ${hisat2_build_memory} GB threshold), so NOT using splice sites and exons to build HISAT2 index.")
        log.info("[HISAT2 index build] Increase memory allocation or lower --hisat2_build_memory to enable splice-aware indexing.")
    }

    """
    mkdir hisat2
    ${extract_exons}
    hisat2-build \\
        -p ${task.cpus} \\
        ${ss} \\
        ${exon} \\
        ${args} \\
        ${fasta} \\
        hisat2/${fasta.baseName}
    """

    stub:
    """
    mkdir hisat2
    """
}
