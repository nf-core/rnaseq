process STRINGTIE_MERGE {
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    path stringtie_gtf
    path annotation_gtf

    output:
    path "stringtie.merged.gtf", emit: gtf
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version'), emit: versions_stringtie, topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    """
    stringtie \\
        --merge $stringtie_gtf \\
        $reference \\
        -o stringtie.merged.gtf \\
        $args

    """

    stub:
    """
    touch stringtie.merged.gtf

    """
}
