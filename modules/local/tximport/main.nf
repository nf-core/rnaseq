process TXIMPORT {
    label "process_medium"

    conda "bioconda::bioconductor-tximeta=1.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.12.0--r41hdfd78af_0' :
        'biocontainers/bioconductor-tximeta:1.12.0--r41hdfd78af_0' }"

    input:
    path ("quants/*")
    path  tx2gene
    val quant_type

    output:
    path "*gene_tpm.tsv"                 , emit: tpm_gene
    path "*gene_counts.tsv"              , emit: counts_gene
    path "*gene_counts_length_scaled.tsv", emit: counts_gene_length_scaled
    path "*gene_counts_scaled.tsv"       , emit: counts_gene_scaled
    path "*gene_lengths.tsv"             , emit: lengths_gene
    path "*transcript_tpm.tsv"           , emit: tpm_transcript
    path "*transcript_counts.tsv"        , emit: counts_transcript
    path "*transcript_lengths.tsv"       , emit: lengths_transcript
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    prefix = task.ext.prefix ?: "${quant_type}.merged"
    """
    tximport.r \\
        NULL \\
        quants \\
        $prefix \\
        $quant_type \\
        $tx2gene

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
