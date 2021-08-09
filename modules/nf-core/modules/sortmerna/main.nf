// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SORTMERNA {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::sortmerna=4.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sortmerna:4.2.0--0"
    } else {
        container "quay.io/biocontainers/sortmerna:4.2.0--0"
    }

    input:
    tuple val(meta), path(reads)
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path   "versions.yml"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def Refs = ""
    for (i=0; i<fasta.size(); i++) { Refs+= " --ref ${fasta[i]}" }
    if (meta.single_end) {
        """
        sortmerna \\
            $Refs \\
            --reads $reads \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --other non_rRNA_reads \\
            $options.args

        gzip -f < non_rRNA_reads.fq > ${prefix}.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            - ${getSoftwareName(task.process)}: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    } else {
        """
        sortmerna \\
            $Refs \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2 \\
            $options.args

        gzip -f < non_rRNA_reads_fwd.fq > ${prefix}_1.fastq.gz
        gzip -f < non_rRNA_reads_rev.fq > ${prefix}_2.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        ${getModuleName(task.process)}:
            - ${getSoftwareName(task.process)}: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    }
}
