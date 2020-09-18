// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path fasta
    path gtf
    val  options

    output:
    path "rsem"         , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    mkdir rsem
    rsem-prepare-reference \\
        --gtf $gtf \\
        --num-threads $task.cpus \\
        $ioptions.args \\
        $fasta \\
        rsem/genome

    rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g" > ${software}.version.txt
    """
}
