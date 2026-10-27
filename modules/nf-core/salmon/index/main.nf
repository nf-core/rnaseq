process SALMON_INDEX {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1ce42a19f9e7135babf14432e80b33ec717a13f14734d150d53347e979919629/data' :
        'community.wave.seqera.io/library/salmon:2.7.0--74784226202c61b9' }"

    input:
    tuple val(meta), path(transcript_fasta), path(genome_fasta)

    output:
    tuple val(meta), path("salmon"), emit: index
    tuple val("${task.process}"), val('salmon'), eval("salmon --version | sed 's/salmon //'"), emit: versions_salmon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def decoys = ''
    def fasta = transcript_fasta
    if (genome_fasta) {
        if ("${genome_fasta}".endsWith('.gz')) {
            genome_fasta = "<(gunzip -c ${genome_fasta})"
        }
        decoys='-d decoys.txt'
        fasta='gentrome.fa'
    }
    if ("${transcript_fasta}".endsWith('.gz')) {
        transcript_fasta = "<(gunzip -c ${transcript_fasta})"
    }
    """
    if [ -n '$genome_fasta' ]; then
        grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 | sed 's/>//g' > decoys.txt
        cat $transcript_fasta $genome_fasta > $fasta
    fi

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $fasta \\
        $decoys \\
        $args \\
        -i salmon
    """

    stub:
    """
    mkdir salmon
    touch salmon/duplicate_clusters.tsv
    touch salmon/index.ctab
    touch salmon/index.ectab
    touch salmon/index.refinfo
    touch salmon/index.ssi
    touch salmon/index.ssi.mphf
    touch salmon/index.tct
    touch salmon/index.tdct
    touch salmon/info.json
    touch salmon/refseq.bin
    touch salmon/refseq_offsets.json
    """
}
