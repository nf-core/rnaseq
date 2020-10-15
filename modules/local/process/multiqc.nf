// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    }

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path fail_mapping_summary
    path fail_strand_check
    path ('fastqc/*')
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
    path ('preseq/*')
    path ('qualimap/*')
    path ('dupradar/*')
    path ('edger/*')
    path ('rseqc/bam_stat/*')
    path ('rseqc/infer_experiment/*')
    path ('rseqc/inner_distance/*')
    path ('rseqc/junction_annotation/*')
    path ('rseqc/junction_saturation/*')
    path ('rseqc/read_distribution/*')
    path ('rseqc/read_duplication/*')
    path ('featurecounts/*')
    path ('featurecounts/biotype/*')
    
    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $options.args $custom_config .
    """
}
