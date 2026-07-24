process TXIMETA_TXIMPORT {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdba33f8ad1b2df156f8f6775279fb217ce0f8233a3dc637337245de9ad2f29f/data' :
        'community.wave.seqera.io/library/bioconductor-tximeta_jq:78bccd386c46a07c' }"

    input:
    tuple val(meta), path("quants/*")
    tuple val(meta2), path(tx2gene)
    val quant_type

    output:
    tuple val(meta), path("*gene_tpm.tsv")                 , emit: tpm_gene
    tuple val(meta), path("*gene_counts.tsv")              , emit: counts_gene
    tuple val(meta), path("*gene_counts_length_scaled.tsv"), emit: counts_gene_length_scaled
    tuple val(meta), path("*gene_counts_scaled.tsv")       , emit: counts_gene_scaled
    tuple val(meta), path("*gene_lengths.tsv")             , emit: lengths_gene
    tuple val(meta), path("*transcript_tpm.tsv")           , emit: tpm_transcript
    tuple val(meta), path("*transcript_counts.tsv")        , emit: counts_transcript
    tuple val(meta), path("*transcript_lengths.tsv")       , emit: lengths_transcript
    tuple val(meta), path("*tx2gene_augmented.tsv")        , emit: tx2gene_augmented
    path "versions.yml"                                    , emit: versions, topic: versions

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
    touch ${meta.id}.tx2gene_augmented.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
