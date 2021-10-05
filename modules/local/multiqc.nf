// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::multiqc=1.10.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.10.1--pyhdfd78af_1"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--pyhdfd78af_1"
    }

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path fail_mapping_summary
    path fail_strand_check
    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')
    path ('sortmerna/*')
    path ('star/*')
    path ('hisat2/*')
    path ('rsem/*')
    path ('salmon/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('picard/markduplicates/*')
    path ('featurecounts/*')
    path ('deseq2/aligner/*')
    path ('deseq2/aligner/*')
    path ('deseq2/pseudoaligner/*')
    path ('deseq2/pseudoaligner/*')
    path ('preseq/*')
    path ('qualimap/*')
    path ('dupradar/*')
    path ('rseqc/bam_stat/*')
    path ('rseqc/infer_experiment/*')
    path ('rseqc/inner_distance/*')
    path ('rseqc/junction_annotation/*')
    path ('rseqc/junction_saturation/*')
    path ('rseqc/read_distribution/*')
    path ('rseqc/read_duplication/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $options.args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
