// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process DESEQ2_QC {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "nfcore/chipseq:1.2.1"
    } else {
        container "nfcore/chipseq:1.2.1"
    }
    
    input:
    path counts
    
    output:    
    path "*.pdf"            , emit: pdf
    path "*.RData"          , emit: rdata
    path "*pca.vals.txt"    , emit: pca_txt
    path "*sample.dists.txt", emit: dists_txt
    path "*log"             , emit: log
    path "size_factors"     , emit: size_factors
    path  "*.version.txt"   , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --cores $task.cpus \\
        $options.args

    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='${software}.version.txt')"
    """
}
