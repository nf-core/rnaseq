// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/rsem:1.3.3--pl526ha52163a_0"
    //container "https://depot.galaxyproject.org/singularity/rsem:1.3.3--pl526ha52163a_0"

    conda (params.conda ? "bioconda::rsem=1.3.3" : null)

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
