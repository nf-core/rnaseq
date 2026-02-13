nextflow.preview.types = true

record DupRadarResult {
    meta:         Map
    scatter:      Path
    boxplot:      Path
    histogram:    Path
    gene_data:    Path
    intercept:    Path
    multiqc:      Path
    session_info: Path
}

process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24bb76357588d05b5637e2954f2dfb3ba04e3eb1ff52c927ffe1906d7d69915a/data' :
        'community.wave.seqera.io/library/bioconductor-dupradar:1.38.0--831da16eb40a64ab' }"

    input:
    (meta: Map, bam: Path): Record
    (meta2: Map, gtf: Path): Record

    output:
    record(
        meta:         meta,
        scatter:      file("*_duprateExpDens.pdf"),
        boxplot:      file("*_duprateExpBoxplot.pdf"),
        histogram:    file("*_expressionHist.pdf"),
        gene_data:    file("*_dupMatrix.txt"),
        intercept:    file("*_intercept_slope.txt"),
        multiqc:      file("*_mqc.txt"),
        session_info: file("*.R_sessionInfo.log")
    )
    tuple val("${task.process}"), val('bioconductor-dupradar'), eval('Rscript -e "library(dupRadar); cat(as.character(packageVersion(\'dupRadar\')))"'), topic: versions

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
    """
}
