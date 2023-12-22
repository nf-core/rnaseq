process HISAT2_BUILD {
    tag "$fasta"
    label 'process_high'
    label 'process_high_memory'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::hisat2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3' :
        'biocontainers/hisat2:2.2.1--h1b792b2_3' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(splicesites)

    output:
    tuple val(meta), path("hisat2") , emit: index
    path "versions.yml"             , emit: versions

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
        extract_exons = gtf ? "hisat2_extract_exons.py $gtf > ${gtf.baseName}.exons.txt" : ""
        ss = splicesites ? "--ss $splicesites" : ""
        exon = gtf ? "--exon ${gtf.baseName}.exons.txt" : ""
    } else {
        log.info "[HISAT2 index build] Less than ${hisat2_build_memory} GB available, so NOT using splice sites and exons to build HISAT2 index."
        log.info "[HISAT2 index build] Use --hisat2_build_memory [small number] to skip this check."
    }
    def VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
