process RSEM_CALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23651ffd6a171ef3ba867cb97ef615f6dd6be39158df9466fe92b5e844cd7d59/data' :
        'community.wave.seqera.io/library/rsem_star:5acb4e8c03239c32' }"

    input:
    tuple val(meta), path(reads)  // FASTQ files or BAM file for --alignments mode
    path  index

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    tuple val(meta), path("*.log")             , emit: logs, optional:true
    path  "versions.yml"                       , emit: versions

    tuple val(meta), path("*.STAR.genome.bam")       , optional:true, emit: bam_star
    tuple val(meta), path("${prefix}.genome.bam")    , optional:true, emit: bam_genome
    tuple val(meta), path("${prefix}.transcript.bam"), optional:true, emit: bam_transcript

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--strandedness forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--strandedness reverse'
    }

    // Detect if input is BAM file(s)
    def is_bam = reads.toString().toLowerCase().endsWith('.bam')
    def alignment_mode = is_bam ? '--alignments' : ''

    // Use metadata for paired-end detection if available, otherwise empty (auto-detect)
    def paired_end = meta.containsKey('single_end') ? (meta.single_end ? "" : "--paired-end") : "unknown"

    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`

    # Use metadata-based paired-end detection, or auto-detect if no metadata provided
    PAIRED_END_FLAG="$paired_end"
    if [ "${paired_end}" == "unknown" ]; then
        # Auto-detect only if no metadata provided
        if [ "${is_bam}" == "true" ]; then
            samtools flagstat $reads | grep -q 'paired in sequencing' && PAIRED_END_FLAG="--paired-end"
        else
            [ ${reads.size()} -gt 1 ] && PAIRED_END_FLAG="--paired-end"
        fi
    fi
    
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --temporary-folder ./tmp/ \\
        $alignment_mode \\
        $strandedness \\
        \$PAIRED_END_FLAG \\
        $args \\
        $reads \\
        \$INDEX \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_bam = reads.toString().toLowerCase().endsWith('.bam')
    """
    touch ${prefix}.genes.results
    touch ${prefix}.isoforms.results
    touch ${prefix}.stat
    touch ${prefix}.log
    
    # Only create STAR BAM output when not in alignment mode
    if [ "${is_bam}" == "false" ]; then
        touch ${prefix}.STAR.genome.bam
    fi
    
    touch ${prefix}.genome.bam
    touch ${prefix}.transcript.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
