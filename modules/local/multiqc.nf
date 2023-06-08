process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path methods_description
    path logo
    path fail_trimming_summary
    path fail_mapping_summary
    path fail_strand_check
    path ('fastqc/raw/*')
    path ('fastqc/trim/*')
    path ('trim_log/*')
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
    path ('rseqc/tin/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
