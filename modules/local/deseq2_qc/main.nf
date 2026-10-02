nextflow.enable.types = true

record Deseq2QcInput {
    id:                        String
    meta:                      Map
    counts_gene_length_scaled: Path
}

record Deseq2Qc {
    id:            String
    meta:          Map
    rdata:         Path?
    pca_vals:      Path?
    plots_pdf:     Path?
    sample_dists:  Path?
    size_factors:  Path?
    log:           Path?
    pca_multiqc:   Path?
    dists_multiqc: Path?
}

process DESEQ2_QC {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d425b12748ce54c44c01a535a1ef5867a6e16cbf62c43151012e893444b1673/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-ggplot2_r-rcolorbrewer_pruned:9e75394d0bc21987' }"

    input:
    sample: Deseq2QcInput
    pca_header_multiqc: Path
    clustering_header_multiqc: Path

    output:
    record(
        id:            sample.id,
        meta:          sample.meta,
        rdata:         file('*.RData',               optional: true),
        pca_vals:      file('*pca.vals.txt',         optional: true),
        plots_pdf:     file('*.pdf',                 optional: true),
        sample_dists:  file('*sample.dists.txt',     optional: true),
        size_factors:  file('size_factors',          optional: true),
        log:           file('*.log',                 optional: true),
        pca_multiqc:   file('*pca.vals_mqc.tsv',     optional: true),
        dists_multiqc: file('*sample.dists_mqc.tsv', optional: true)
    ) as Deseq2Qc

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'deseq2_qc.r'

    stub:
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    prefix = task.ext.prefix ?: "deseq2"
    """
    touch ${label_lower}.pca.vals_mqc.tsv
    touch ${label_lower}.sample.dists_mqc.tsv
    touch ${prefix}.dds.RData
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.plots.pdf
    touch ${prefix}.sample.dists.txt
    touch R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS

    mkdir size_factors
    touch size_factors/${prefix}.size_factors.RData
    # One per-sample size_factors file per data column in ${sample.counts_gene_length_scaled}; the
    # module test snaps these names so the stub must mirror real-run output.
    for i in `head ${sample.counts_gene_length_scaled} -n 1 | cut -f3-`;
    do
        touch size_factors/\${i}.size_factors.RData
    done
    """
}
