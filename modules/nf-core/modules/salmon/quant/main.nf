// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SALMON_QUANT {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::salmon=1.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.4.0--hf69c8f4_0"
    } else {
        container "quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    path  transcript_fasta
    val   alignment_mode
    val   lib_type

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path   "versions.yml"             , emit: version

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def reference   = "--index $index"
    def input_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    if (alignment_mode) {
        reference   = "-t $transcript_fasta"
        input_reads = "-a $reads"
    }

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $options.args \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getModuleName(task.process)}:
        $software: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
