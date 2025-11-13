process BBMAP_BBSPLIT {
    tag "$meta.id"
    label 'process_high'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  primary_ref
    tuple val(other_ref_names), path(other_ref_paths)
    val   only_build_index

    output:
    path "bbsplit"                            , optional:true, emit: index
    tuple val(meta), path('*primary*fastq.gz'), optional:true, emit: primary_fastq
    tuple val(meta), path('*fastq.gz')        , optional:true, emit: all_fastq
    tuple val(meta), path('*txt')             , optional:true, emit: stats
    tuple val(meta), path('*.log')            , optional:true, emit: log
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

    def fastq_in=''
    def fastq_out=''
    def index_files=''
    def refstats_cmd=''

    if (only_build_index) {
        if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = 'ref_primary=' +primary_ref + ' ' + other_refs.join(' ') + ' path=bbsplit'
        } else {
            log.error 'ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files.'
        }
    } else {
        if (index) {
            index_files = "path=$index"
        } else if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = "ref_primary=${primary_ref} ${other_refs.join(' ')}"
        } else {
            log.error 'ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files.'
        }
        fastq_in  = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
        fastq_out = meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"
        refstats_cmd = 'refstats=' + prefix + '.stats.txt'
    }
    """

    # When we stage in the index files the time stamps get disturbed, which
    # bbsplit doesn't like. Fix the time stamps in its summaries. This needs to
    # be done via Java to match what bbmap does

    if [ $index ]; then
        for summary_file in \$(find $index/ref/genome -name summary.txt); do
            src=\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/bbsplit|bbsplit|')
            mod=\$(echo "System.out.println(java.nio.file.Files.getLastModifiedTime(java.nio.file.Paths.get(\\"\$src\\")).toMillis());" | jshell -J-Djdk.lang.Process.launchMechanism=vfork -)
            sed "s|^last modified.*|last modified\\t\$mod|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
        done
    fi

    # Run BBSplit

    bbsplit.sh \\
        -Xmx${avail_mem}M \\
        $index_files \\
        threads=$task.cpus \\
        $fastq_in \\
        $fastq_out \\
        $refstats_cmd \\
        $args 2>| >(tee ${prefix}.log >&2)

    # Summary files will have an absolute path that will make the index
    # impossible to use in other processes- we can fix that

    for summary_file in \$(find bbsplit/ref/genome -name summary.txt); do
        src=\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/bbsplit|bbsplit|')
        sed "s|^source.*|source\\t\$src|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def other_refs = ''
    other_ref_names.eachWithIndex { name, index ->
        other_refs += "echo '' | gzip > ${prefix}_${name}.fastq.gz"
    }
    """
    if [ ! -d bbsplit ]; then
        mkdir bbsplit
    fi

    if ! (${only_build_index}); then
        echo '' | gzip >  ${prefix}_primary.fastq.gz
        ${other_refs}
        touch ${prefix}.stats.txt
    fi

    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
