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
    path  index, name: 'input_index'
    path  primary_ref
    tuple val(other_ref_names), path(other_ref_paths)
    val   only_build_index

    output:
    path "bbsplit_index"                      , optional:true, emit: index
    tuple val(meta), path('*primary*fastq.gz'), optional:true, emit: primary_fastq
    tuple val(meta), path('*fastq.gz')        , optional:true, emit: all_fastq
    tuple val(meta), path('*txt')             , optional:true, emit: stats
    tuple val(meta), path('*.log')            , optional:true, emit: log
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), topic: versions, emit: versions_bbmap

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
    other_ref_names.eachWithIndex { name, idx ->
        other_refs << "ref_${name}=${other_ref_paths[idx]}"
    }

    def fastq_in=''
    def fastq_out=''
    def index_files=''
    def refstats_cmd=''
    def use_index = index ? true : false

    if (only_build_index) {
        if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = 'ref_primary=' +primary_ref + ' ' + other_refs.join(' ') + ' path=bbsplit_build'
        } else {
            log.error 'ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files.'
        }
    } else {
        if (index) {
            index_files = "path=index_writable"
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

    # If using a pre-built index, create writable structure: symlink all files except
    # summary.txt (which we copy to modify). When we stage in the index files the time
    # stamps get disturbed, which bbsplit doesn't like. Fix the time stamps in summaries.
    if [ "$use_index" == "true" ]; then
        find input_index/ref -type f | while read -r f; do
            target="index_writable/\${f#input_index/}"
            mkdir -p "\$(dirname "\$target")"
            [[ \$(basename "\$f") == "summary.txt" ]] && cp "\$f" "\$target" || ln -s "\$(realpath "\$f")" "\$target"
        done
        find index_writable/ref/genome -name summary.txt | while read -r summary_file; do
            src=\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/ref/|index_writable/ref/|')
            mod=\$(echo "System.out.println(java.nio.file.Files.getLastModifiedTime(java.nio.file.Paths.get(\\"\$src\\")).toMillis());" | jshell -J-Djdk.lang.Process.launchMechanism=vfork - 2>/dev/null | grep -oE '^[0-9]{12,14}\$')
            sed -e 's|bbsplit_index/ref|index_writable/ref|' -e "s|^last modified.*|last modified\\t\$mod|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
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
    # impossible to use in other processes - fix paths and rename atomically
    if [ -d bbsplit_build/ref/genome ]; then
        find bbsplit_build/ref/genome -name summary.txt | while read -r summary_file; do
            sed "s|^source.*|source\\t\$(grep '^source' "\$summary_file" | cut -f2- -d\$'\\t' | sed 's|.*/bbsplit_build|bbsplit_index|')|" "\$summary_file" > \${summary_file}.tmp && mv \${summary_file}.tmp \${summary_file}
        done
        mv bbsplit_build bbsplit_index
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def other_refs = ''
    other_ref_names.eachWithIndex { name, _idx ->
        other_refs += "echo '' | gzip > ${prefix}_${name}.fastq.gz"
    }
    def will_build_index = only_build_index || (!index && primary_ref && other_ref_names && other_ref_paths)
    """
    # Create index directory if building an index (either only_build_index or on-the-fly)
    if [ "${will_build_index}" == "true" ]; then
        mkdir -p bbsplit_index
    fi

    # Only create output files if splitting (not just building index)
    if ! (${only_build_index}); then
        echo '' | gzip >  ${prefix}_primary.fastq.gz
        ${other_refs}
        touch ${prefix}.stats.txt
    fi

    touch ${prefix}.log
    """
}
