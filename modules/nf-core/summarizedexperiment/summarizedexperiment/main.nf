process SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    meta            : Map
    matrix_files    : List<Path>
    rowdata         : Path
    coldata         : Path

    output:
    rds : Path = file("*.rds")
    log : Path = file("*.R_sessionInfo.log")

    script:
    template 'summarizedexperiment.r'

    stub:
    """
    touch ${meta.id}.SummarizedExperiment.rds
    touch ${meta.id}.R_sessionInfo.log
    """
}
