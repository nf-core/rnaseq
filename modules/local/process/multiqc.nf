// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    //container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"

    conda (params.conda ? "bioconda::multiqc=1.9" : null)

    input:
    path multiqc_config
    path mqc_custom_config
    path software_versions
    path workflow_summary
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
    val options

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "multiqc_plots"       , emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def ioptions      = initOptions(options)
    def rtitle        = custom_runName ? "--title \"$custom_runName\"" : ''
    def rfilename     = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def custom_config = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $ioptions.args $rtitle $rfilename $custom_config .
    """
}
