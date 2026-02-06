process SENTIEON_RSEMCALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39a3e1a85912520836ad054c8ac0497b463bb5170e0e907183dbd08509dad997/data' :
        'community.wave.seqera.io/library/rsem_sentieon:3e4315fa0b636313' }"

    input:
    tuple val(meta), path(reads)  // FASTQ files or BAM file for --alignments mode
    path  index

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    tuple val(meta), path("*.log")             , emit: logs, optional:true

    tuple val(meta), path("*.STAR.genome.bam")       , optional:true, emit: bam_star
    tuple val(meta), path("${prefix}.genome.bam")    , optional:true, emit: bam_genome
    tuple val(meta), path("${prefix}.transcript.bam"), optional:true, emit: bam_transcript

    tuple val("${task.process}"), val('rsem'), eval('rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g"'), topic: versions, emit: versions_rsem
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions, emit: versions_star
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version 2>&1 | sed -e "s/sentieon-genomics-//g"'), topic: versions, emit: versions_sentieon

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

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""

    """
    $sentieonLicense

    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`

    # Create symlink to sentieon in PATH
    ln -sf \$(which sentieon) ./STAR
    export PATH=".:\$PATH"

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
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_bam = reads.toString().toLowerCase().endsWith('.bam')
    """
    # Create symlink to sentieon in PATH for version detection
    ln -sf \$(which sentieon) ./STAR
    export PATH=".:\$PATH"

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
    """
}
