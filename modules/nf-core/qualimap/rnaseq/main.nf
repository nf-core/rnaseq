process QUALIMAP_RNASEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.3--hdfd78af_0' :
        'biocontainers/qualimap:2.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('qualimap'), eval("qualimap 2>&1 | sed -n 's/.*QualiMap v.\\(.*\\)/\\1/p'"), emit: versions_qualimap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    def memory = (task.memory.mega*0.8).intValue() + 'M'

    def strandedness = 'non-strand-specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }
    """
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        rnaseq \\
        $args \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        $paired_end \\
        -outdir $prefix
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    """
}
