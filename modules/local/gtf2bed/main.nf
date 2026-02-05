process GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    path gtf

    output:
    path '*.bed', emit: bed
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | sed -n 's/.*v\\(.*\\)) built.*/\\1/p'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    gtf2bed \\
        $gtf \\
        > ${gtf.baseName}.bed
    """

    stub:
    """
    touch ${gtf.baseName}.bed
    """
}
