nextflow.enable.types = true

process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
        : 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'}"

    input:
    record(id: String, meta: Map, reads: List<Path>)

    stage:
    stageAs reads, '?/*'

    output:
    record(
        id:   id,
        meta: meta,
        html: files('*.html').toSorted { f -> f.name },
        zip:  files('*.zip').toSorted { f -> f.name }
    )

    topic:
    tuple(task.process, 'fastqc', eval('fastqc --version | sed "/FastQC v/!d; s/.*v//"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // New names for the symlinks made in the bash while loop below
    def new_names = reads.size() == 1
        ? ["${prefix}.${reads[0].extension}"]
        : reads.withIndex().collect { entry, index -> "${prefix}_${index + 1}.${entry.extension}" }.toList()
    def rename_to = reads.withIndex().collect { entry, index -> "${entry} ${new_names[index]}" }.join(' ')
    def renamed_files = new_names.join(' ')

    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpus allows to stick to requested amount of RAM in the label
    def memory_in_mb = task.memory
        ? (task.memory.toUnit('MB') / task.cpus).intValue()
        : null
    // FastQC memory value allowed range (100 - 10000)
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    def fastqc_memory_arg = fastqc_memory ? "--memory ${fastqc_memory}" : ''

    """
    printf "%s %s\\n" ${rename_to} | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${fastqc_memory_arg} \\
        ${renamed_files}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip
    """
}
