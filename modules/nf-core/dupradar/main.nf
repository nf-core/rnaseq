nextflow.enable.types = true

process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24bb76357588d05b5637e2954f2dfb3ba04e3eb1ff52c927ffe1906d7d69915a/data' :
        'community.wave.seqera.io/library/bioconductor-dupradar:1.38.0--831da16eb40a64ab' }"

    input:
    tuple(meta: Map, bam: Path)
    tuple(meta2: Map, gtf: Path)

    output:
    record(
        meta:            meta,
        scatter2d:       file("*_duprateExpDens.pdf"),
        boxplot:         file("*_duprateExpBoxplot.pdf"),
        hist:            file("*_expressionHist.pdf"),
        dupmatrix:       file("*_dupMatrix.txt"),
        intercept_slope: file("*_intercept_slope.txt"),
        multiqc:         files("*_mqc.txt", optional: true).toSorted { f -> f.name },
        session_info:    file("*.R_sessionInfo.log")
    )

    topic:
    file('versions.yml') >> 'versions'

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
