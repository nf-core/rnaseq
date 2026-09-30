nextflow.enable.types = true

record SalmonIndexInput {
    id:               String
    meta:             Map
    transcript_fasta: Path
    genome_fasta:     Path?
}

record SalmonIndexResult {
    id:    String
    meta:  Map
    index: Path
}

process SALMON_INDEX {
    tag "$sample.meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1ce42a19f9e7135babf14432e80b33ec717a13f14734d150d53347e979919629/data' :
        'community.wave.seqera.io/library/salmon:2.7.0--74784226202c61b9' }"

    input:
    sample: SalmonIndexInput

    output:
    record(id: sample.id, meta: sample.meta, index: file('salmon'))

    topic:
    tuple(task.process, 'salmon', eval("salmon --version | sed 's/salmon //'")) >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def decoys = ''
    def fasta = "${sample.transcript_fasta}"
    def genome = sample.genome_fasta ? "${sample.genome_fasta}" : ''
    def transcripts = "${sample.transcript_fasta}"
    if (sample.genome_fasta) {
        if ("${sample.genome_fasta}".endsWith('.gz')) {
            genome = "<(gunzip -c ${sample.genome_fasta})"
        }
        decoys='-d decoys.txt'
        fasta='gentrome.fa'
    }
    if ("${sample.transcript_fasta}".endsWith('.gz')) {
        transcripts = "<(gunzip -c ${sample.transcript_fasta})"
    }
    """
    if [ -n '$genome' ]; then
        grep '^>' $genome | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 | sed 's/>//g' > decoys.txt
        cat $transcripts $genome > $fasta
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
