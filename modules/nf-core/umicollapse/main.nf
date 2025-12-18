process UMICOLLAPSE {
    tag "$meta.id"
    label "process_high"
    label "process_high_memory"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umicollapse:1.1.0--hdfd78af_0' :
        'biocontainers/umicollapse:1.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(bai)
    val(mode)

    output:
    tuple val(meta), path("*.bam"),                   emit: bam,        optional: true
    tuple val(meta), path("*dedup*fastq.gz"),         emit: fastq,      optional: true
    tuple val(meta), path("*_UMICollapse.log"),       emit: log
    path  "versions.yml" ,                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.0-0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // Memory allocation: We need to make sure that both heap and stack size is sufficiently large for
    // umicollapse. We set the stack size to 5% of the available memory, the heap size to 90%
    // which leaves 5% for stuff happening outside of java without the scheduler killing the process.
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999 //most java jdks will not allow Xss > 1GB, so fixing this to the allowed max
    if ( mode !in [ 'fastq', 'bam' ] ) {
        error "Mode must be one of 'fastq' or 'bam'."
    }
    extension = mode.contains("fastq") ? "fastq.gz" : "bam"
    """
    # Getting the umicollapse jar file like this because `umicollapse` is a Python wrapper script generated
    # by conda that allows to set the heap size (Xmx), but not the stack size (Xss).
    # `which` allows us to get the directory that contains `umicollapse`, independent of whether we
    # are in a container or conda environment.
    UMICOLLAPSE_JAR=\$(dirname \$(which umicollapse))/../share/umicollapse-${VERSION}/umicollapse.jar
    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M \\
        -jar \$UMICOLLAPSE_JAR \\
        $mode \\
        -i ${input} \\
        -o ${prefix}.${extension} \\
        $args | tee ${prefix}_UMICollapse.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0-1'
    if ( mode !in [ 'fastq', 'bam' ] ) {
        error "Mode must be one of 'fastq' or 'bam'."
    }
    extension = mode.contains("fastq") ? "fastq.gz" : "bam"
    """
    touch ${prefix}.dedup.${extension}
    touch ${prefix}_UMICollapse.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: $VERSION
    END_VERSIONS
    """
}
