def VERSION = '2.2.0' // Version information not provided by tool on CLI

process HISAT2_BUILD {
    tag "$fasta"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? 'bioconda::hisat2=2.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3' :
        'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3' }"

    input:
    path fasta
    path gtf
    path splicesites

    output:
    path "hisat2"       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 0
    if (!task.memory) {
        log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
    } else {
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        avail_mem = task.memory.toGiga()
    }

    def ss = ''
    def exon = ''
    def extract_exons = ''
    def hisat2_build_memory = params.hisat2_build_memory ? (params.hisat2_build_memory as nextflow.util.MemoryUnit).toGiga() : 0
    if (avail_mem >= hisat2_build_memory) {
        log.info "[HISAT2 index build] At least ${hisat2_build_memory} GB available, so using splice sites and exons to build HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.exons.txt"
        ss = "--ss $splicesites"
        exon = "--exon ${gtf.baseName}.exons.txt"
    } else {
        log.info "[HISAT2 index build] Less than ${hisat2_build_memory} GB available, so NOT using splice sites and exons to build HISAT2 index."
        log.info "[HISAT2 index build] Use --hisat2_build_memory [small number] to skip this check."
    }
    """
    mkdir hisat2
    $extract_exons
    hisat2-build \\
        -p $task.cpus \\
        $ss \\
        $exon \\
        $args \\
        $fasta \\
        hisat2/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: $VERSION
    END_VERSIONS
    """
}
