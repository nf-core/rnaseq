// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEM_CALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rsem:1.3.3--pl526ha52163a_0"
    //container "https://depot.galaxyproject.org/singularity/rsem:1.3.3--pl526ha52163a_0"

    conda (params.conda ? "bioconda::rsem=1.3.3" : null)

    input:
    tuple val(meta), path(reads)
    path  index
    val   options

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    path  "*.version.txt"                      , emit: version

    tuple val(meta), path("*.STAR.genome.bam")       , optional:true, emit: bam_star
    tuple val(meta), path("${prefix}.genome.bam")    , optional:true, emit: bam_genome
    tuple val(meta), path("${prefix}.transcript.bam"), optional:true, emit: bam_transcript

    script:
    def software   = getSoftwareName(task.process)
    def ioptions   = initOptions(options)
    prefix         = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def paired_end = meta.single_end ? "" : "--paired-end"
    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/.grp//'`
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --temporary-folder ./tmp/
        $paired_end \\
        $ioptions.args \\
        $reads \\
        $INDEX \\
        $prefix

    rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g" > ${software}.version.txt
    """
}
