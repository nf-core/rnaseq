process BBMAP_BBSPLIT {
    tag "$meta.id"
    label 'process_high'
    label 'error_retry'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  primary_ref
    tuple val(other_ref_names), path (other_ref_paths)
    val   only_build_index

    output:
    path "bbsplit"                            , optional:true, emit: index
    tuple val(meta), path('*primary*fastq.gz'), optional:true, emit: primary_fastq
    tuple val(meta), path('*fastq.gz')        , optional:true, emit: all_fastq
    tuple val(meta), path('*txt')             , optional:true, emit: stats
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[BBSplit] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def other_refs = []
    other_ref_names.eachWithIndex { name, index ->
        other_refs << "ref_${name}=${other_ref_paths[index]}"
    }
    if (only_build_index) {
        if (primary_ref && other_ref_names && other_ref_paths) {
            """
            bbsplit.sh \\
                -Xmx${avail_mem}M \\
                ref_primary=$primary_ref \\
                ${other_refs.join(' ')} \\
                path=bbsplit \\
                threads=$task.cpus \\
                $args

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
            END_VERSIONS
            """
        } else {
            log.error 'ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files.'
        }
    } else {
        def index_files = ''
        if (index) {
            index_files = "path=$index"
        } else if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = "ref_primary=${primary_ref} ${other_refs.join(' ')}"
        } else {
            log.error 'ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files.'
        }
        def fastq_in  = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
        def fastq_out = meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"
        """
        bbsplit.sh \\
            -Xmx${avail_mem}M \\
            $index_files \\
            threads=$task.cpus \\
            $fastq_in \\
            $fastq_out \\
            refstats=${prefix}.stats.txt \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        END_VERSIONS
        """
    }
}
