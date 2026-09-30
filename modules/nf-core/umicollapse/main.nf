nextflow.enable.types = true

process UMICOLLAPSE {
    tag "${meta.id}"
    label "process_high"
    label "process_high_memory"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/umicollapse:1.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/umicollapse:1.1.0--hdfd78af_0'}"

    input:
    tuple(meta: Map, input: Path, bai: Path?)
    mode: String

    output:
    record(
        id:    meta.id,
        meta:  meta,
        bam:   file('*.bam', optional: true),
        fastq: file('*dedup*fastq.gz', optional: true),
        log:   file('*_UMICollapse.log')
    )

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple(task.process, 'umicollapse', '1.1.0-0') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Memory allocation: We need to make sure that both heap and stack size is sufficiently large for
    // umicollapse. We set the stack size to 5% of the available memory, the heap size to 90%
    // which leaves 5% for stuff happening outside of java without the scheduler killing the process.
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999
    //most java jdks will not allow Xss > 1GB, so fixing this to the allowed max
    if (mode !in ['fastq', 'bam']) {
        error("Mode must be one of 'fastq' or 'bam'.")
    }
    def extension = mode.contains("fastq") ? "fastq.gz" : "bam"
    """
    # The generated launcher allows configuring heap size, but not stack size.
    UMICOLLAPSE_JAR=\$(find "\$(dirname "\$(command -v umicollapse)")/../share" -maxdepth 2 -name umicollapse.jar -print -quit)

    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M \\
        -jar "\$UMICOLLAPSE_JAR" \\
        ${mode} \\
        -i ${input} \\
        -o ${prefix}.${extension} \\
        ${args} | tee ${prefix}_UMICollapse.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (mode !in ['fastq', 'bam']) {
        error("Mode must be one of 'fastq' or 'bam'.")
    }
    def extension = mode.contains("fastq") ? "fastq.gz" : "bam"
    """
    touch ${prefix}.dedup.${extension}
    touch ${prefix}_UMICollapse.log
    """
}
