nextflow.enable.types = true

record SummarizedexperimentInput {
    id:           String
    meta:         Map
    matrix_files: List<Path>
}

record SummarizedexperimentSummarizedexperimentResult {
    id:   String
    meta: Map
    rds:  Path
    log:  Path
}

process SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    sample: SummarizedexperimentInput
    rowdata: Path?
    coldata: Path?

    output:
    record(
        id:   sample.id,
        meta: sample.meta,
        rds:  file("*.rds"),
        log:  file("*.R_sessionInfo.log")
    ) as SummarizedexperimentSummarizedexperimentResult

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'summarizedexperiment.r'

    stub:
    """
    touch ${sample.meta.id}.SummarizedExperiment.rds
    touch ${sample.meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
