nextflow.enable.types = true

process SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    tuple(meta: Map, matrix_files: List<Path>)
    tuple(meta2: Map, rowdata: Path?)
    tuple(meta3: Map, coldata: Path?)

    output:
    record(
        id:   meta.id,
        meta: meta,
        rds:  file("*.rds"),
        log:  file("*.R_sessionInfo.log")
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'summarizedexperiment.r'

    stub:
    """
    touch ${meta.id}.SummarizedExperiment.rds
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
