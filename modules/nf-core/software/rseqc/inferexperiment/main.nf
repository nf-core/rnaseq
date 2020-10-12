// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    //container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"

    conda (params.conda ? "bioconda::rseqc=3.0.1" : null)

    input:
    tuple val(meta), path(bam)
    path  bed
    val   options

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt
    path  "*.version.txt"                          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        $ioptions.args \\
        > ${prefix}.infer_experiment.txt

    infer_experiment.py --version | sed -e "s/infer_experiment.py //g" > ${software}.version.txt
    """
}
