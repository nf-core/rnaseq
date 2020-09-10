// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    // Don't upgrade me - 2.7X indices incompatible with iGenomes.
    container "quay.io/biocontainers/star:2.6.1d--0"
    //container "https://depot.galaxyproject.org/singularity/star:2.6.1d--0"

    conda (params.conda ? "bioconda::star=2.6.1d" : null)

    input:
    path fasta
    path gtf
    val  options

    output:
    path "star"         , emit: index
    path "*.version.txt", emit: version

    script:
    def software  = getSoftwareName(task.process)
    def ioptions  = initOptions(options)
    def memory    = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --runThreadN $task.cpus \\
        $memory \\
        $ioptions.args

    STAR --version | sed -e "s/STAR_//g" > ${software}.version.txt
    """
}
