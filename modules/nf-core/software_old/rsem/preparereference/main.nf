// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::rsem=1.3.3 bioconda::star=2.7.6a" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    }

    input:
    path fasta
    path gtf
    
    output:
    path "rsem"         , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    mkdir rsem
    rsem-prepare-reference \\
        --gtf $gtf \\
        --num-threads $task.cpus \\
        $options.args \\
        $fasta \\
        rsem/genome

    rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g" > ${software}.version.txt
    """
}
