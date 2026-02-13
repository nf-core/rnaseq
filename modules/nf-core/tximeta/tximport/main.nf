nextflow.preview.types = true

record TximportResult {
    meta:                       Map
    tpm_gene:                   Path
    counts_gene:                Path
    counts_gene_length_scaled:  Path
    counts_gene_scaled:         Path
    lengths_gene:               Path
    tpm_transcript:             Path
    counts_transcript:          Path
    lengths_transcript:         Path
}

process TXIMETA_TXIMPORT {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta%3A1.20.1--r43hdfd78af_0' :
        'biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_0' }"

    input:
    (meta: Map, quants: Path): Record
    (meta2: Map, tx2gene: Path): Record

    stage:
    stageAs(quants, 'quants/*')
    quant_type: String

    output:
    record(
        meta:                      meta,
        tpm_gene:                  file("*gene_tpm.tsv"),
        counts_gene:               file("*gene_counts.tsv"),
        counts_gene_length_scaled: file("*gene_counts_length_scaled.tsv"),
        counts_gene_scaled:        file("*gene_counts_scaled.tsv"),
        lengths_gene:              file("*gene_lengths.tsv"),
        tpm_transcript:            file("*transcript_tpm.tsv"),
        counts_transcript:         file("*transcript_counts.tsv"),
        lengths_transcript:        file("*transcript_lengths.tsv")
    )
    tuple val("${task.process}"), val('bioconductor-tximeta'), eval('Rscript -e "library(tximeta); cat(as.character(packageVersion(\'tximeta\')))"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tximport.r'

    stub:
    """
    touch ${meta.id}.gene_tpm.tsv
    touch ${meta.id}.gene_counts.tsv
    touch ${meta.id}.gene_counts_length_scaled.tsv
    touch ${meta.id}.gene_counts_scaled.tsv
    touch ${meta.id}.gene_lengths.tsv
    touch ${meta.id}.transcript_tpm.tsv
    touch ${meta.id}.transcript_counts.tsv
    touch ${meta.id}.transcript_lengths.tsv
    """
}
