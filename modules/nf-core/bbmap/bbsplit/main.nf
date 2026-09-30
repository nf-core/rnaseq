nextflow.enable.types = true

include { ReadsInput } from '../../types'

process BBMAP_BBSPLIT {
    tag "${sample.meta.id}"
    label 'process_high'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    sample: ReadsInput
    index: Path?
    primary_ref: Path?
    tuple(other_ref_names: List<String>, other_ref_paths: List<Path>)
    only_build_index: Boolean

    stage:
    stageAs index, 'input_index'

    output:
    record(
        id:                 sample.id,
        meta:               sample.meta,
        index:              file('bbsplit_index', optional: true),
        reads:              files('*primary*fastq.gz', optional: true).toSorted { f -> f.name },
        other_genome_reads: files('*fastq.gz', optional: true).findAll { f -> !f.name.contains('primary') }.toSorted { f -> f.name },
        stats:              file('*txt', optional: true),
        log:                file('*.log', optional: true)
    )

    topic:
    tuple(task.process, 'bbmap', eval('bbversion.sh | grep -v "Duplicate cpuset"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[BBSplit] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.toMega()*0.8).intValue()
    }

    def other_refs = other_ref_names.withIndex().collect { name, idx -> "ref_${name}=${other_ref_paths[idx]}" }

    def fastq_in=''
    def fastq_out=''
    def index_files=''
    def refstats_cmd=''
    def use_index = index ? true : false

    if (only_build_index) {
        if (primary_ref && other_ref_names && other_ref_paths) {
            index_files = "ref_primary=${primary_ref} ${other_refs.join(' ')} path=bbsplit_build"
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
        fastq_in  = sample.meta.single_end ? "in=${sample.reads[0]}" : "in=${sample.reads[0]} in2=${sample.reads[1]}"
        fastq_out = sample.meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"
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
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def other_refs = other_ref_names.collect { name -> "echo '' | gzip > ${prefix}_${name}.fastq.gz" }.join('')
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
