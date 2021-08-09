// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RSEQC_INNERDISTANCE {
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
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript
    path   "versions.yml"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            $options.args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt

        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            - ${getSoftwareName(task.process)}: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    } else {
        """
        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            - ${getSoftwareName(task.process)}: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    }
}
