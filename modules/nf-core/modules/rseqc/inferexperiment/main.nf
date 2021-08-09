// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::rseqc=3.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt
    path   "versions.yml"                          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        $options.args \\
        > ${prefix}.infer_experiment.txt

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """
}
