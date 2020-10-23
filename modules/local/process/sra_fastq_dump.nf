// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Download SRA data via parallel-fastq-dump
 */
process SRA_FASTQ_DUMP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::parallel-fastq-dump=0.6.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/parallel-fastq-dump:0.6.6--py_1"
    } else {
        container "quay.io/biocontainers/sra-tools:2.10.8--pl526haddd2b5_0"
        // container "quay.io/biocontainers/parallel-fastq-dump:0.6.6--py_1"
        //container "quay.io/biocontainers/parallel-fastq-dump:0.6.6--py_0"
        // container "quay.io/biocontainers/parallel-fastq-dump:0.6.5--py_0"
        // container "nfcore/viralrecon:1.1.0"
    }
    
    input:
    val meta

    // output:
    // tuple val(meta), path("*fastq.gz"), emit: fastq
    // tuple val(meta), path("*log")     , emit: log
    // path "*.version.txt"              , emit: version

    script:
    id         = "${meta.id.split('_')[0..-2].join('_')}"
    paired_end = meta.single_end ? "" : "--readids --split-e"
    rm_orphan  = meta.single_end ? "" : "[ -f  ${id}.fastq.gz ] && rm ${id}.fastq.gz"    
    """
    prefetch --progress $id
    vdb-validate --verbose $id
    """
}
// mkdir -p ~/.ncbi
// printf '/LIBS/GUID = "%s"\n' `uuid` > ~/.ncbi/user-settings.mkfg
//export HOME=./
//vdb-config --interactive

    // parallel-fastq-dump \\
    //     --sra-id $id \\
    //     --threads $task.cpus \\
    //     --outdir ./ \\
    //     --tmpdir ./ \\
    //     --gzip \\
    //     $paired_end \\
    //     > ${id}.fastq_dump.log
    // $rm_orphan

    // echo \$(parallel-fastq-dump --version 2>&1) | sed 's/^.*parallel-fastq-dump : //; s/ "fastq-dump".*\$//' > parallel-fastq-dump.version.txt
