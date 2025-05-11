process TXIMETA_TXIMPORT {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta%3A1.20.1--r43hdfd78af_0' :
        'biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_0' }"

    input:
    meta        : Map
    quants      : List<Path>
    tx2gene     : Path
    quant_type  : String

    stage:
    stageAs "quants/*", quants

    output:
    tpm_gene                    : Path = file("*gene_tpm.tsv")
    counts_gene                 : Path = file("*gene_counts.tsv")
    counts_gene_length_scaled   : Path = file("*gene_counts_length_scaled.tsv")
    counts_gene_scaled          : Path = file("*gene_counts_scaled.tsv")
    lengths_gene                : Path = file("*gene_lengths.tsv")
    tpm_transcript              : Path = file("*transcript_tpm.tsv")
    counts_transcript           : Path = file("*transcript_counts.tsv")
    lengths_transcript          : Path = file("*transcript_lengths.tsv")

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
