// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getModuleName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? "bioconda::star=2.6.1d bioconda::samtools=1.10 conda-forge::gawk=5.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0"
    }

    input:
    path fasta
    path gtf

    output:
    path "star"         , emit: index
    path  "versions.yml", emit: version

    script:
    def software = getSoftwareName(task.process)
    def memory   = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def args     = options.args.tokenize()
    if (args.contains('--genomeSAindexNbases')) {
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $options.args

        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            $software: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        samtools faidx $fasta
        NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            $options.args

        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            $software: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }
}
