process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24bb76357588d05b5637e2954f2dfb3ba04e3eb1ff52c927ffe1906d7d69915a/data' :
        'community.wave.seqera.io/library/bioconductor-dupradar:1.38.0--831da16eb40a64ab' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*_duprateExpDens.pdf")   , emit: scatter2d
    tuple val(meta), path("*_duprateExpBoxplot.pdf"), emit: boxplot
    tuple val(meta), path("*_expressionHist.pdf")   , emit: hist
    tuple val(meta), path("*_dupMatrix.txt")        , emit: dupmatrix
    tuple val(meta), path("*_intercept_slope.txt")  , emit: intercept_slope
    tuple val(meta), path("*_mqc.txt")              , emit: multiqc
    tuple val(meta), path("*.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dupradar.r'

    stub:
    """
    touch ${meta.id}_duprateExpDens.pdf
    touch ${meta.id}_duprateExpBoxplot.pdf
    touch ${meta.id}_expressionHist.pdf
    touch ${meta.id}_dupMatrix.txt
    touch ${meta.id}_intercept_slope.txt
    touch ${meta.id}_dup_intercept_mqc.txt
    touch ${meta.id}_duprateExpDensCurve_mqc.txt
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-dupradar: \$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
    END_VERSIONS
    """
}
