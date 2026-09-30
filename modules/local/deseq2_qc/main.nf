nextflow.enable.types = true

process DESEQ2_QC {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d425b12748ce54c44c01a535a1ef5867a6e16cbf62c43151012e893444b1673/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-ggplot2_r-rcolorbrewer_pruned:9e75394d0bc21987' }"

    input:
    record(id: String, meta: Map, counts_gene_length_scaled: Path)
    pca_header_multiqc: Path
    clustering_header_multiqc: Path

    output:
    record(
        id:            id,
        meta:          meta,
        rdata:         file('*.RData',               optional: true),
        pca_vals:      file('*pca.vals.txt',         optional: true),
        plots_pdf:     file('*.pdf',                 optional: true),
        sample_dists:  file('*sample.dists.txt',     optional: true),
        size_factors:  file('size_factors',          optional: true),
        log:           file('*.log',                 optional: true),
        pca_multiqc:   file('*pca.vals_mqc.tsv',     optional: true),
        dists_multiqc: file('*sample.dists_mqc.tsv', optional: true)
    )

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'
    tuple(task.process, 'bioconductor-deseq2', eval("Rscript -e \"library(DESeq2); cat(as.character(packageVersion('DESeq2')))\"")) >> 'versions'

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()
    prefix = task.ext.prefix ?: "deseq2"
    """
    deseq2_qc.r \\
        --count_file $counts_gene_length_scaled \\
        --outdir ./ \\
        --cores $task.cpus \\
        --outprefix $prefix \\
        $args

    if [ -f "R_sessionInfo.log" ]; then
        # Handle PCA files
        sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" <$pca_header_multiqc > pca_header.tmp
        sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" pca_header.tmp
        cat pca_header.tmp *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv
        rm pca_header.tmp

        # Handle clustering files
        sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" <$clustering_header_multiqc > clustering_header.tmp
        sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" clustering_header.tmp
        cat clustering_header.tmp *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
        rm clustering_header.tmp
    fi
    """

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

    mkdir size_factors
    touch size_factors/${prefix}.size_factors.RData
    # One per-sample size_factors file per data column in $counts_gene_length_scaled; the
    # module test snaps these names so the stub must mirror real-run output.
    for i in `head $counts_gene_length_scaled -n 1 | cut -f3-`;
    do
        touch size_factors/\${i}.size_factors.RData
    done
    """
}
