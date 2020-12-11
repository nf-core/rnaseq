// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
params.multiqc_label = ''
def options          = initOptions(params.options)

process DESEQ2_QC {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.2 conda-forge::r-rcolorbrewer=1.1_2 conda-forge::r-pheatmap=1.0.12 bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel=1.22.0 bioconda::bioconductor-tximport=1.16.0 bioconda::bioconductor-complexheatmap=2.4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
    }

    input:
    path counts
    path pca_header_multiqc
    path clustering_header_multiqc
    
    output:    
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path  "*.version.txt"       , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def label_lower = params.multiqc_label.toLowerCase()
    def label_upper = params.multiqc_label.toUpperCase()
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --cores $task.cpus \\
        $options.args
    
    if [ -f "R_sessionInfo.log" ]; then
        sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" <$pca_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
        cat tmp.txt *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv

        sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" <$clustering_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
        cat tmp.txt *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
    fi

    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='${software}.version.txt')"
    """
}