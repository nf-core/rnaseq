// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0"
    //container "https://depot.galaxyproject.org/singularity/rseqc:4.0.0--py38h0213d0e_0"

    conda (params.conda ? "bioconda::rseqc=4.0.0" : null)

    input:
    tuple val(meta), path(bam)
    path  bed
    val   options

    output:
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript
    path  "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    if (!meta.single_end) {
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            $ioptions.args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt

        inner_distance.py --version | sed -e "s/inner_distance.py //g" > ${software}.version.txt
        """
    } else {
        """
        inner_distance.py --version | sed -e "s/inner_distance.py //g" > ${software}.version.txt
        """
    }
}
