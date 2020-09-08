// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process SALMON_QUANT {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/salmon:1.3.0--hf69c8f4_0"
    //container "https://depot.galaxyproject.org/singularity/salmon:1.3.0--hf69c8f4_0"

    conda (params.conda ? "bioconda::salmon=1.3.0" : null)

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    val   options

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    prefix       = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def strandedness = meta.single_end ? 'U' : 'IU'
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? 'SF' : 'ISF'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? 'SR' : 'ISR'
    }
    def endedness = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        --index $index \\
        $endedness \\
        $ioptions.args \\
        -o $prefix

    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}
