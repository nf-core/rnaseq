process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.32.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-dupradar:1.32.0--r43hdfd78af_0' }"

    input:
    meta    : Map
    bam     : Path
    gtf     : Path

    output:
    scatter2d       : Path = file("*_duprateExpDens.pdf")
    boxplot         : Path = file("*_duprateExpBoxplot.pdf")
    hist            : Path = file("*_expressionHist.pdf")
    dupmatrix       : Path = file("*_dupMatrix.txt")
    intercept_slope : Path = file("*_intercept_slope.txt")
    session_info    : Path = file("*.R_sessionInfo.log")

    topic:
    file('versions.yml') >> 'versions'
    file('*_mqc.txt') >> 'logs'

    script:
    template 'dupradar.r'

    stub:
    """
    touch ${meta.id}_duprateExpDens.pdf
    touch ${meta.id}_duprateExpBoxplot.pdf
    touch ${meta.id}_expressionHist.pdf
    touch ${meta.id}_dupMatrix.txt
    touch ${meta.id}_intercept_slope.txt
    touch ${meta.id}_mqc.txt
    touch ${meta.id}.R_sessionInfo.log
    """
}
