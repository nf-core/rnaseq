process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2' :
        'biocontainers/salmon:1.10.3--h6dccd9a_2' }"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon"      , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    if (genome_fasta.endsWith('.gz')) {
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
        gentrome      = "gentrome.fa.gz"
    }
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decoys.txt \\
        $args \\
        -i salmon

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir salmon
    touch salmon/complete_ref_lens.bin
    touch salmon/ctable.bin
    touch salmon/ctg_offsets.bin
    touch salmon/duplicate_clusters.tsv
    touch salmon/info.json
    touch salmon/mphf.bin
    touch salmon/pos.bin
    touch salmon/pre_indexing.log
    touch salmon/rank.bin
    touch salmon/refAccumLengths.bin
    touch salmon/ref_indexing.log
    touch salmon/reflengths.bin
    touch salmon/refseq.bin
    touch salmon/seq.bin
    touch salmon/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
