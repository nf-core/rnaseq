// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BBMAP_BBSPLIT {
    //tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bbmap=38.93" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:38.93--he522d1c_0"
    } else {
        container "quay.io/biocontainers/bbmap:38.93--he522d1c_0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  fastas
    val   index_only

    output:
    path "bbsplit"                    , optional:true, emit: index
    tuple val(meta), path('*fastq.gz'), optional:true, emit: fastq
    tuple val(meta), path('*txt')     , optional:true, emit: stats
    path "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[BBSplit] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    def fasta_list = []
    fastas.eachWithIndex { fasta, index ->
        fasta_list << "ref_${index}=${fasta}"
    }
    if (index_only) {
        """
        bbsplit.sh \\
            -Xmx${avail_mem}g \\
            ${fasta_list.join(' ')} \\
            path=bbsplit \\
            threads=$task.cpus \\
            $options.args

        echo \$(bbversion.sh) > ${software}.version.txt
        """
    } else {
        def fastq_in  = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
        def fastq_out = meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"
        """
        bbsplit.sh \\
            -Xmx${avail_mem}g \\
            ${fasta_list.join(' ')} \\
            path=$index \\
            threads=$task.cpus \\
            $fastq_in \\
            $fastq_out \\
            refstats=${prefix}.stats.txt \\
            $options.args

        echo \$(bbversion.sh) > ${software}.version.txt
        """
    }
}
