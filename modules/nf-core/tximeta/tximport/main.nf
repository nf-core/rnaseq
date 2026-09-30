nextflow.enable.types = true

include { QuantsInput } from '../../types'

process TXIMETA_TXIMPORT {
    tag "${sample.meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdba33f8ad1b2df156f8f6775279fb217ce0f8233a3dc637337245de9ad2f29f/data' :
        'community.wave.seqera.io/library/bioconductor-tximeta_jq:78bccd386c46a07c' }"

    input:
    sample: QuantsInput
    tx2gene: Path
    quant_type: String

    stage:
    stageAs sample.quants, 'quants/*'

    output:
    record(
        id:                        sample.id,
        meta:                      sample.meta,
        tpm_gene:                  file("*gene_tpm.tsv"),
        counts_gene:               file("*gene_counts.tsv"),
        lengths_gene:              file("*gene_lengths.tsv"),
        counts_gene_length_scaled: file("*gene_counts_length_scaled.tsv"),
        counts_gene_scaled:        file("*gene_counts_scaled.tsv"),
        tpm_transcript:            file("*transcript_tpm.tsv"),
        counts_transcript:         file("*transcript_counts.tsv"),
        lengths_transcript:        file("*transcript_lengths.tsv"),
        tx2gene_augmented:         file("*tx2gene_augmented.tsv")
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'tximport.r'

    stub:
    """
    touch ${sample.meta.id}.gene_tpm.tsv
    touch ${sample.meta.id}.gene_counts.tsv
    touch ${sample.meta.id}.gene_counts_length_scaled.tsv
    touch ${sample.meta.id}.gene_counts_scaled.tsv
    touch ${sample.meta.id}.gene_lengths.tsv
    touch ${sample.meta.id}.transcript_tpm.tsv
    touch ${sample.meta.id}.transcript_counts.tsv
    touch ${sample.meta.id}.transcript_lengths.tsv
    touch ${sample.meta.id}.tx2gene_augmented.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
