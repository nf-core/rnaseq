process GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'biocontainers/perl:5.26.2' }"

    input:
    gtf : Path

    output:
    bed : Path = file('*.bed')

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    gtf2bed \\
        $gtf \\
        > ${gtf.baseName}.bed
    """
}
