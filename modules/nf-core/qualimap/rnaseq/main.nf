nextflow.enable.types = true

process QUALIMAP_RNASEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.3--hdfd78af_0' :
        'quay.io/biocontainers/qualimap:2.3--hdfd78af_0' }"

    input:
    tuple(meta: Map, bam: Path)
    tuple(meta2: Map, gtf: Path)

    output:
    record(meta: meta, qualimap: file("${prefix}"))

    topic:
    tuple(task.process, 'qualimap', eval("qualimap -h | sed -n 's/^QualiMap v.//p'")) >> 'versions'

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    def memory = "${(task.memory.toMega() * 0.8).intValue()}M"

    """
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=${memory} \\
        rnaseq \\
        ${args} \\
        -bam ${bam} \\
        -gtf ${gtf} \\
        ${paired_end} \\
        -outdir ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    """
}
